"""
Kamodo Vector Field Tracer - Core field line tracing functionality

High-performance vector field tracer for Kamodo models.
Supports magnetic field lines, velocity streamlines, current streamlines, etc.
All coordinates are in Earth radii (RE).

Optimized with:
- Numba JIT compilation for core numerical routines
- Parallel processing using joblib for multiple field lines
- Vectorized operations for improved performance
- Memory optimization techniques
- Batch processing for efficient handling of large datasets
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import warnings
from functools import lru_cache
import gc

# Import numba components with error handling
try:
    from numba import jit, prange
    from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
    NUMBA_AVAILABLE = True
    # Suppress numba warnings about deprecated NumPy API features
    warnings.filterwarnings("ignore", category=NumbaDeprecationWarning)
    warnings.filterwarnings("ignore", category=NumbaPendingDeprecationWarning)
except ImportError:
    print("Warning: Numba not available. Performance optimizations will be limited.")
    NUMBA_AVAILABLE = False
    # Create dummy jit decorator
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator
    prange = range

# Import joblib with error handling
try:
    from joblib import Parallel, delayed
    JOBLIB_AVAILABLE = True
except ImportError:
    print("Warning: Joblib not available. Parallel processing will be disabled.")
    JOBLIB_AVAILABLE = False

# Import tqdm if available for progress bars
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    # Create a simple fallback tqdm function
    def tqdm(iterable=None, **kwargs):
        return iterable
    TQDM_AVAILABLE = False


def get_colormap(name):
    """
    Get colormap in a way that is compatible with different matplotlib versions.

    Parameters:
    -----------
    name : str
        Name of the colormap

    Returns:
    --------
    cmap : matplotlib colormap
        The requested colormap
    """
    try:
        # Try new interface first (matplotlib >= 3.7)
        import matplotlib
        if hasattr(matplotlib, 'colormaps'):
            return matplotlib.colormaps[name]
        else:
            # Fallback for matplotlib < 3.7
            return plt.cm.get_cmap(name)
    except (AttributeError, KeyError):
        # Additional fallback
        try:
            return plt.get_cmap(name)
        except:
            # Last resort - return a default colormap
            print(f"Warning: Could not find colormap '{name}', using 'viridis'")
            return plt.get_cmap('viridis')


@jit(nopython=True, fastmath=True, cache=True)
def _rk4_step_numba(y, h, f_vec):
    """
    Perform a single 4th-order Runge-Kutta step, optimized with numba.
    
    Parameters:
    -----------
    y : array
        Current position
    h : float
        Step size
    f_vec : array
        Field vector at current position
    
    Returns:
    --------
    y_new : array
        New position after step
    """
    # For a full RK4 implementation, we'd need multiple field evaluations
    # This is a simplified version for use with external field calls
    k1 = h * f_vec
    return y + k1


@jit(nopython=True, fastmath=True, cache=True)
def _compute_length_numba(coordinates):
    """
    Compute total length of field line path using numba.
    
    Parameters:
    -----------
    coordinates : array
        Array of shape (n_points, 3) containing trace coordinates
    
    Returns:
    --------
    length : float
        Total length of field line in RE
    """
    length = 0.0
    for i in range(1, len(coordinates)):
        dx = coordinates[i, 0] - coordinates[i-1, 0]
        dy = coordinates[i, 1] - coordinates[i-1, 1]
        dz = coordinates[i, 2] - coordinates[i-1, 2]
        length += np.sqrt(dx*dx + dy*dy + dz*dz)
    return length


@jit(nopython=True, fastmath=True, cache=True)
def _normalize_vector_numba(vector):
    """
    Normalize a vector using numba optimization.
    
    Parameters:
    -----------
    vector : array
        Input vector
        
    Returns:
    --------
    normalized : array
        Normalized vector
    magnitude : float
        Original magnitude
    """
    magnitude = np.sqrt(np.sum(vector * vector))
    if magnitude > 0:
        return vector / magnitude, magnitude
    else:
        return vector, 0.0


@jit(nopython=True, fastmath=True, cache=True)
def _check_bounds_numba(position, min_r, max_r):
    """
    Check if position is within bounds using numba.
    
    Parameters:
    -----------
    position : array
        Position vector [x, y, z]
    min_r : float
        Minimum allowed radius
    max_r : float
        Maximum allowed radius
        
    Returns:
    --------
    within_bounds : bool
        True if within bounds
    """
    r = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
    return min_r <= r <= max_r


class KamodoVectorFieldTracer:
    """
    Class to trace vector field lines in Kamodo models.
    
    Attributes:
    -----------
    kamodo_object : object
        Kamodo object with vector field components
    field_type : str
        Type of field to trace ('magnetic', 'velocity', 'current', etc.)
    component_names : list
        List of field component variable names in kamodo object
    use_numba : bool
        Whether to use numba optimization
    use_parallel : bool
        Whether to use parallel processing
    n_jobs : int
        Number of parallel jobs for multi-line tracing
    """
    
    def __init__(self, kamodo_object, field_type='magnetic',
                 component_names=None, use_numba=True, use_parallel=False, n_jobs=-1,
                 verbose=False):
        """
        Initialize vector field tracer.
        
        Parameters:
        -----------
        kamodo_object : object
            Kamodo object with vector field components
        field_type : str
            Type of field to trace ('magnetic', 'velocity', 'current', etc.)
        component_names : list, optional
            List of field component variable names in kamodo object.
            If None, default names based on field_type will be used.
        use_numba : bool, optional (default=True)
            Whether to use numba optimization for core routines
        use_parallel : bool, optional (default=False)
            Whether to use parallel processing for multiple lines
        n_jobs : int, optional (default=-1)
            Number of jobs for parallel processing (-1 uses all cores)
        verbose : bool, optional (default=False)
            Whether to print detailed information during execution
        """
        self.kamodo_object = kamodo_object
        self.field_type = field_type.lower()
        self.use_numba = use_numba and NUMBA_AVAILABLE
        self.use_parallel = use_parallel and JOBLIB_AVAILABLE
        self.n_jobs = n_jobs
        self.verbose = verbose
        
        # Set up component names based on field type
        if component_names is None:
            if self.field_type == 'magnetic':
                prefix = 'B'
            elif self.field_type == 'velocity':
                prefix = 'v'
            elif self.field_type == 'current':
                prefix = 'J'
            elif self.field_type == 'electric':
                prefix = 'E'
            else:
                prefix = 'f'  # Generic field
            
            self.component_names = [f'{prefix}_x', f'{prefix}_y', f'{prefix}_z']
        else:
            self.component_names = component_names
        
        # Cache field component functions for faster access
        self._cache_field_components()
        
        # Default tracing parameters
        self.default_params = {
            'max_steps': 5000,
            'step_size_re': 0.05,
            'direction': 'both',
            'min_altitude_re': 1.0,
            'max_altitude_re': 50.0,
            'adaptive_stepping': True,
            'batch_size': 100,  # For batch processing
            'max_memory_mb': 2000  # For memory management
        }
        
        # Performance tracking
        self.performance_stats = {
            'total_traces': 0,
            'total_time': 0.0,
            'total_points': 0
        }
        
        # For Kamodo calling convention
        self._determine_call_method()
        
        if self.verbose:
            self._print_configuration()
    
    def _print_configuration(self):
        """Print the current configuration of the tracer."""
        print(f"\nKamodoVectorFieldTracer Configuration:")
        print(f"  Field type: {self.field_type}")
        print(f"  Component names: {self.component_names}")
        print(f"  Numba optimization: {'enabled' if self.use_numba else 'disabled'}")
        print(f"  Parallel processing: {'enabled' if self.use_parallel else 'disabled'}")
        print(f"  Number of jobs: {self.n_jobs}")
        print(f"  Default parameters: {self.default_params}")
        print(f"  Kamodo calling method: {self._call_method}")
    
    def _cache_field_components(self):
        """Cache field component functions for faster access."""
        self.field_components = []
        
        available_vars = list(self.kamodo_object.keys()) if hasattr(self.kamodo_object, 'keys') else []
        
        for component in self.component_names:
            if component in self.kamodo_object:
                self.field_components.append(self.kamodo_object[component])
            else:
                raise ValueError(f"Component '{component}' not found in Kamodo object. "
                                 f"Available variables: {available_vars}")
    
    def _determine_call_method(self):
        """
        Determine the correct method to call Kamodo functions.
        Kamodo functions typically accept [time, x, y, z] arrays.
        """
        # Default calling method (works with most Kamodo models)
        self._call_method = 'array_with_time'
        
        if self.verbose:
            print(f"Using Kamodo calling method: {self._call_method}")
            print(f"For format: component([time, x, y, z])")
    
    def get_field_vector(self, position_re, time_hours):
        """
        Get field vector at given position and time.
        Uses Kamodo calling convention: component([time, x, y, z])
        
        Parameters:
        -----------
        position_re : array_like
            Position in RE [x, y, z]
        time_hours : float
            Time in hours
        
        Returns:
        --------
        field_vector : ndarray
            Field vector at given position [Bx, By, Bz] or [vx, vy, vz] etc.
        field_magnitude : float
            Magnitude of field vector
        """
        x, y, z = position_re
        
        # Kamodo calling format: component([time, x, y, z])
        coords = [time_hours, x, y, z]
        
        # Get field components
        field_components_values = []
        for i, component in enumerate(self.field_components):
            try:
                value = component(coords)
                field_components_values.append(float(value))
            except Exception as e:
                raise ValueError(f"Error evaluating {self.component_names[i]} at "
                               f"coords [t={time_hours:.3f}, x={x:.3f}, y={y:.3f}, z={z:.3f}]: {e}")
        
        # Convert to numpy array
        field_vector = np.array(field_components_values, dtype=np.float64)
        
        # Calculate field magnitude
        if self.use_numba:
            _, field_magnitude = _normalize_vector_numba(field_vector)
        else:
            field_magnitude = np.sqrt(np.sum(field_vector**2))
        
        return field_vector, field_magnitude
    
    def get_field_vectors_batch(self, positions_re, time_hours):
        """
        Get field vectors for multiple positions efficiently using Kamodo's batch capability.
        
        Parameters:
        -----------
        positions_re : array_like
            Array of shape (n_points, 3) containing positions in RE
        time_hours : float
            Time in hours
        
        Returns:
        --------
        field_vectors : ndarray
            Array of shape (n_points, 3) containing field vectors
        field_magnitudes : ndarray
            Array of shape (n_points,) containing field magnitudes
        """
        positions_re = np.array(positions_re)
        if positions_re.ndim == 1:
            positions_re = positions_re.reshape(1, -1)
        
        n_points = len(positions_re)
        
        # Create coordinate array for Kamodo: [[time, x, y, z], [time, x, y, z], ...]
        coords_array = []
        for i in range(n_points):
            x, y, z = positions_re[i]
            coords_array.append([time_hours, x, y, z])
        
        # Get field components for all points
        field_components_values = []
        for i, component in enumerate(self.field_components):
            try:
                # Handle both single and multiple points
                if n_points == 1:
                    values = [component(coords_array[0])]
                else:
                    values = component(coords_array)
                
                # Ensure we have a list/array of values
                if not isinstance(values, (list, np.ndarray)):
                    values = [values]
                
                field_components_values.append(np.array(values, dtype=np.float64))
                
            except Exception as e:
                raise ValueError(f"Error evaluating {self.component_names[i]} for batch of {n_points} points: {e}")
        
        # Stack into (n_points, 3) array - transpose to get [Bx,By,Bz] for each point
        field_vectors = np.column_stack(field_components_values)
        
        # Calculate field magnitudes
        if self.use_numba:
            field_magnitudes = np.array([_normalize_vector_numba(vec)[1] for vec in field_vectors])
        else:
            field_magnitudes = np.sqrt(np.sum(field_vectors**2, axis=1))
        
        return field_vectors, field_magnitudes
    
    def _check_stopping_criteria(self, position_re, steps_taken, field_mag):
        """
        Check if field line tracing should stop.
        
        Parameters:
        -----------
        position_re : array_like
            Current position in RE [x, y, z]
        steps_taken : int
            Number of steps taken
        field_mag : float
            Field magnitude
        
        Returns:
        --------
        stop : bool
            Whether to stop tracing
        reason : str
            Reason for stopping
        """
        # Check bounds using optimized function if available
        if self.use_numba:
            within_bounds = _check_bounds_numba(
                np.array(position_re), 
                self.params['min_altitude_re'], 
                self.params['max_altitude_re']
            )
            if not within_bounds:
                r = np.sqrt(np.sum(np.array(position_re)**2))
                if r < self.params['min_altitude_re']:
                    return True, "hit_min_altitude"
                else:
                    return True, "hit_max_altitude"
        else:
            # Calculate distance from origin (Earth center)
            r = np.linalg.norm(position_re)
            
            # Check stopping criteria
            if r < self.params['min_altitude_re']:
                return True, "hit_min_altitude"
            elif r > self.params['max_altitude_re']:
                return True, "hit_max_altitude"
        
        # Other stopping criteria
        if steps_taken >= self.params['max_steps']:
            return True, "max_steps_reached"
        elif field_mag < 1e-10:  # Field too weak
            return True, "field_too_weak"
        
        return False, "none"

    def trace_vector_line(self, start_position_re, time_hours, **kwargs):
        """
        Trace a vector field line starting from given position.
        
        Parameters:
        -----------
        start_position_re : array_like
            Starting position in RE [x, y, z]
        time_hours : float
            Time in hours
        **kwargs : dict
            Additional parameters:
                - max_steps: Maximum number of steps
                - step_size_re: Initial step size in RE
                - direction: 'forward', 'backward', or 'both'
                - min_altitude_re: Minimum altitude in RE
                - max_altitude_re: Maximum altitude in RE
                - adaptive_stepping: Whether to use adaptive step sizing
                - batch_size: Size for batch processing
                - max_memory_mb: Maximum memory usage in MB
        
        Returns:
        --------
        results : dict
            Dictionary with traced field line results
        """
        # Set parameters from kwargs with defaults
        self.params = self.default_params.copy()
        self.params.update(kwargs)
        
        start_time = time.time()
        
        # Convert start position to numpy array
        start_position_re = np.array(start_position_re, dtype=np.float64)
        
        # Determine which directions to trace
        trace_forward = self.params['direction'] in ['forward', 'both']
        trace_backward = self.params['direction'] in ['backward', 'both']
        
        results = {}
        
        # Trace in forward direction
        if trace_forward:
            forward = self._trace_single_direction(
                start_position_re, time_hours, forward=True)
            results['forward'] = forward
        
        # Trace in backward direction
        if trace_backward:
            backward = self._trace_single_direction(
                start_position_re, time_hours, forward=False)
            results['backward'] = backward
        
        # Combine results for bidirectional tracing
        if trace_forward and trace_backward:
            # Reverse backward coordinates and concatenate
            rev_backward_coords = backward['coordinates'][::-1]
            combined_coords = np.vstack((rev_backward_coords[:-1], forward['coordinates']))
            
            # Combine field magnitudes and components
            rev_backward_field_mag = backward['field_magnitude'][::-1]
            combined_field_mag = np.concatenate((rev_backward_field_mag[:-1], 
                                                forward['field_magnitude']))
            
            rev_backward_field_comp = backward['field_components'][::-1]
            combined_field_comp = np.vstack((rev_backward_field_comp[:-1], 
                                           forward['field_components']))
            
            # Calculate total length
            combined_length = backward['length'] + forward['length']
            
            # Store combined results
            results['combined'] = {
                'coordinates': combined_coords,  # In RE
                'field_magnitude': combined_field_mag,
                'field_components': combined_field_comp,
                'length': combined_length,  # In RE
                'n_points': len(combined_coords),
                'stop_reasons': [backward['stop_reason'], forward['stop_reason']]
            }
        
        # Calculate execution time
        execution_time = time.time() - start_time
        results['execution_time'] = execution_time
        
        # Update performance stats
        self.performance_stats['total_traces'] += 1
        self.performance_stats['total_time'] += execution_time
        if 'combined' in results:
            self.performance_stats['total_points'] += results['combined']['n_points']
        elif 'forward' in results:
            self.performance_stats['total_points'] += results['forward']['n_points']
        
        return results

    def _trace_single_direction(self, start_position_re, time_hours, forward=True):
        """
        Trace field line in single direction from start position.
        Uses correct Kamodo calling convention.
        
        Parameters:
        -----------
        start_position_re : array_like
            Starting position in RE [x, y, z]
        time_hours : float
            Time in hours
        forward : bool
            Whether to trace in forward direction
        
        Returns:
        --------
        result : dict
            Dictionary with traced field line results in single direction
        """
        # Initialize storage
        coordinates = [start_position_re.copy()]
        field_magnitudes = []
        field_components = []
        
        # Get initial field vector
        field_vec, field_mag = self.get_field_vector(start_position_re, time_hours)
        field_magnitudes.append(field_mag)
        field_components.append(field_vec.copy())
        
        # Set direction sign
        sign = 1.0 if forward else -1.0
        
        # Current position
        position = start_position_re.copy()
        
        # Step size (can be adaptive)
        step_size = self.params['step_size_re']
        
        stop = False
        stop_reason = "none"
        steps_taken = 0
        
        # Main integration loop
        while not stop:
            # Normalize field vector to get unit direction
            if field_mag > 0:
                if self.use_numba:
                    unit_vec, _ = _normalize_vector_numba(field_vec)
                else:
                    unit_vec = field_vec / field_mag
            else:
                stop = True
                stop_reason = "zero_field"
                break
            
            # Use adaptive step size if enabled
            if self.params['adaptive_stepping']:
                # Adjust step size based on field magnitude and position
                r = np.linalg.norm(position)
                # Smaller steps near Earth, larger steps far away
                step_size = max(0.01, min(0.2, 0.05 * (r/3.0)))
                
                # Smaller steps in regions of high field curvature
                if len(coordinates) > 2:
                    prev_unit_vec = (field_components[-2] / field_magnitudes[-2] 
                                   if field_magnitudes[-2] > 0 else unit_vec)
                    angle = np.arccos(np.clip(np.dot(unit_vec, prev_unit_vec), -1.0, 1.0))
                    if angle > 0.1:
                        step_size *= 0.5  # Reduce step size in high curvature regions
            
            # Calculate step using RK4 method
            try:
                if self.use_numba:
                    # Simplified Euler step for numba optimization
                    step_vec = sign * step_size * unit_vec
                    new_position = position + step_vec
                else:
                    # Full RK4 integration
                    k1 = sign * step_size * unit_vec
                    
                    # k2 calculation
                    pos_k1_half = position + 0.5 * k1
                    field_vec_k1, field_mag_k1 = self.get_field_vector(pos_k1_half, time_hours)
                    if field_mag_k1 > 0:
                        k2 = sign * step_size * (field_vec_k1 / field_mag_k1)
                    else:
                        k2 = k1
                    
                    # k3 calculation
                    pos_k2_half = position + 0.5 * k2
                    field_vec_k2, field_mag_k2 = self.get_field_vector(pos_k2_half, time_hours)
                    if field_mag_k2 > 0:
                        k3 = sign * step_size * (field_vec_k2 / field_mag_k2)
                    else:
                        k3 = k2
                    
                    # k4 calculation
                    pos_k3 = position + k3
                    field_vec_k3, field_mag_k3 = self.get_field_vector(pos_k3, time_hours)
                    if field_mag_k3 > 0:
                        k4 = sign * step_size * (field_vec_k3 / field_mag_k3)
                    else:
                        k4 = k3
                    
                    # Combined step
                    new_position = position + (k1 + 2*k2 + 2*k3 + k4) / 6.0
            
            except Exception as e:
                print(f"Error during integration step: {e}")
                stop = True
                stop_reason = "integration_error"
                break
            
            # Update position
            position = new_position
            
            # Get field at new position
            try:
                field_vec, field_mag = self.get_field_vector(position, time_hours)
            except Exception as e:
                print(f"Error evaluating field at new position: {e}")
                stop = True
                stop_reason = "field_evaluation_error"
                break
            
            # Store results
            coordinates.append(position.copy())
            field_magnitudes.append(field_mag)
            field_components.append(field_vec.copy())
            
            # Increment step counter
            steps_taken += 1
            
            # Check stopping criteria
            stop, stop_reason = self._check_stopping_criteria(position, steps_taken, field_mag)
        
        # Convert lists to arrays
        coordinates = np.array(coordinates)
        field_magnitudes = np.array(field_magnitudes)
        field_components = np.array(field_components)
        
        # Calculate path length
        if self.use_numba and len(coordinates) > 1:
            length = _compute_length_numba(coordinates)
        else:
            length = 0.0
            for i in range(1, len(coordinates)):
                length += np.linalg.norm(coordinates[i] - coordinates[i-1])
        
        # Return results
        return {
            'coordinates': coordinates,
            'field_magnitude': field_magnitudes,
            'field_components': field_components,
            'length': length,
            'n_points': len(coordinates),
            'stop_reason': stop_reason
        }

    def trace_multiple_lines(self, start_points_re, time_hours, show_progress=False, use_batches=None, **kwargs):
        """
        Trace multiple vector field lines/streamlines.
        Optimized with parallel processing and batch handling for large datasets.

        Parameters:
        -----------
        start_points_re : array_like
            Array of shape (n_lines, 3) containing starting points in RE
        time_hours : float
            Time in hours since midnight of first day of data
        show_progress : bool, optional (default=False)
            Whether to show progress bar
        use_batches : bool, optional
            Whether to process in batches for large datasets (auto-determined if None)
        **kwargs : dict
            Additional arguments passed to trace_vector_line

        Returns:
        --------
        traces : list
            List of trace result dictionaries
        """
        start_points_re = np.array(start_points_re)
        n_lines = len(start_points_re)
        
        # Set parameters from kwargs with defaults
        params = self.default_params.copy()
        params.update(kwargs)
        
        # Determine if we should use batch processing
        if use_batches is None:
            # Auto-determine based on number of points
            use_batches = n_lines > params['batch_size']
        
        if self.verbose:
            print(f"Tracing {n_lines} {self.field_type} field lines at t={time_hours:.2f} hours...")
            print(f"Batch processing: {'enabled' if use_batches else 'disabled'}")
            print(f"Parallel processing: {'enabled' if self.use_parallel else 'disabled'}")
            print(f"Numba optimization: {'enabled' if self.use_numba else 'disabled'}")
        
        start_time = time.time()
        
        # Define tracing function for a single point
        def trace_single(idx, point):
            try:
                return self.trace_vector_line(point, time_hours, **kwargs)
            except Exception as e:
                print(f"Error tracing line {idx}: {e}")
                return None
        
        # Use batch processing for large datasets
        if use_batches:
            traces = self._trace_in_batches(
                start_points_re, 
                trace_single, 
                batch_size=params['batch_size'],
                show_progress=show_progress,
                max_memory_mb=params['max_memory_mb']
            )
        else:
            # Standard processing (parallel or sequential)
            traces = self._trace_without_batches(
                start_points_re,
                trace_single,
                show_progress=show_progress
            )
        
        # Record total execution time
        total_time = time.time() - start_time
        if self.verbose:
            print(f"Completed {len(traces)} field lines in {total_time:.2f} seconds "
                  f"({len(traces) / total_time:.2f} lines/sec)")
        
        return traces
    
    def _trace_in_batches(self, start_points_re, trace_func, batch_size=100, 
                         show_progress=False, max_memory_mb=2000):
        """
        Process field line tracing in batches to manage memory usage.
        
        Parameters:
        -----------
        start_points_re : array_like
            Array of start points
        trace_func : callable
            Function to trace single line
        batch_size : int
            Size of each batch
        show_progress : bool
            Whether to show progress bars
        max_memory_mb : float
            Maximum memory usage target in MB
            
        Returns:
        --------
        traces : list
            Combined list of trace results
        """
        n_points = len(start_points_re)
        n_batches = (n_points + batch_size - 1) // batch_size
        
        if self.verbose:
            print(f"Processing {n_points} points in {n_batches} batches of size {batch_size}")
            print(f"Target max memory usage: {max_memory_mb} MB")
        
        all_results = []
        
        # Create outer progress bar for batches
        batch_iter = range(n_batches)
        if show_progress and TQDM_AVAILABLE:
            batch_iter = tqdm(batch_iter, desc="Processing batches", unit="batch")
        
        for i in batch_iter:
            start_idx = i * batch_size
            end_idx = min((i + 1) * batch_size, n_points)
            batch = start_points_re[start_idx:end_idx]
            
            if self.verbose:
                print(f"Processing batch {i+1}/{n_batches} with {len(batch)} points")
            
            # Process this batch (either in parallel or sequentially)
            batch_results = self._trace_without_batches(
                batch, 
                trace_func, 
                show_progress=show_progress and not TQDM_AVAILABLE  # Only show inner progress if no outer tqdm
            )
            
            all_results.extend(batch_results)
            
            # Memory management - force garbage collection every few batches
            if i % 3 == 0 and i > 0:
                gc.collect()
                if self.verbose:
                    print(f"  Performed garbage collection after batch {i+1}")
        
        return all_results
    
    def _trace_without_batches(self, start_points_re, trace_func, show_progress=False):
        """
        Process field line tracing without batching (either parallel or sequential).
        
        Parameters:
        -----------
        start_points_re : array_like
            Array of start points
        trace_func : callable
            Function to trace single line
        show_progress : bool
            Whether to show progress bar
            
        Returns:
        --------
        traces : list
            List of trace results
        """
        n_lines = len(start_points_re)
        
        # Use parallel processing if enabled and multiple lines
        if self.use_parallel and n_lines > 1:
            # Create a wrapper function for parallel processing
            parallel_func = delayed(trace_func)
            
            # Execute in parallel
            if show_progress and TQDM_AVAILABLE:
                # For parallel processing with tqdm, we need a different approach
                print(f"Running parallel processing on {n_lines} lines...")
                traces = Parallel(n_jobs=self.n_jobs, prefer="threads")(
                    parallel_func(i, point) for i, point in 
                    tqdm(enumerate(start_points_re), total=n_lines, desc="Tracing lines")
                )
            else:
                traces = Parallel(n_jobs=self.n_jobs, prefer="threads")(
                    parallel_func(i, point) for i, point in enumerate(start_points_re)
                )
        else:
            # Sequential processing
            traces = []
            iterator = enumerate(start_points_re)
            
            if show_progress and TQDM_AVAILABLE:
                iterator = tqdm(iterator, total=n_lines, desc="Tracing lines")
                
            for i, point in iterator:
                trace_result = trace_func(i, point)
                if trace_result is not None:
                    traces.append(trace_result)
        
        # Filter out None results (from errors)
        traces = [t for t in traces if t is not None]
        
        return traces
    
    def get_performance_stats(self):
        """
        Get performance statistics for this tracer instance.
        
        Returns:
        --------
        stats : dict
            Dictionary containing performance statistics
        """
        stats = self.performance_stats.copy()
        if stats['total_traces'] > 0:
            stats['avg_time_per_trace'] = stats['total_time'] / stats['total_traces']
            stats['avg_points_per_trace'] = stats['total_points'] / stats['total_traces']
            stats['points_per_second'] = stats['total_points'] / stats['total_time'] if stats['total_time'] > 0 else 0
            stats['traces_per_second'] = stats['total_traces'] / stats['total_time'] if stats['total_time'] > 0 else 0
        else:
            stats['avg_time_per_trace'] = 0
            stats['avg_points_per_trace'] = 0
            stats['points_per_second'] = 0
            stats['traces_per_second'] = 0
        
        return stats
    
    def benchmark(self, start_points_re, time_hours, n_repeats=3, **kwargs):
        """
        Comprehensive benchmark of tracer performance with different optimization settings.
        
        Parameters:
        -----------
        start_points_re : array_like
            Array of shape (n_lines, 3) containing starting points in RE
        time_hours : float
            Time in hours
        n_repeats : int
            Number of times to repeat each benchmark for reliability
        **kwargs : dict
            Additional arguments passed to trace_multiple_lines
        
        Returns:
        --------
        results : dict
            Dictionary containing benchmark results
        """
        print(f"Starting benchmark with {len(start_points_re)} start points, {n_repeats} repeats each")
        
        results = {}
        
        # Save original settings
        orig_numba = self.use_numba
        orig_parallel = self.use_parallel
        
        # Test configurations
        configs = [
            {"name": "baseline", "use_numba": False, "use_parallel": False},
        ]
        
        if NUMBA_AVAILABLE:
            configs.append({"name": "numba_only", "use_numba": True, "use_parallel": False})
            
        if JOBLIB_AVAILABLE:
            configs.append({"name": "parallel_only", "use_numba": False, "use_parallel": True})
            
        if NUMBA_AVAILABLE and JOBLIB_AVAILABLE:
            configs.append({"name": "numba_and_parallel", "use_numba": True, "use_parallel": True})
        
        for config in configs:
            # Set configuration
            self.use_numba = config["use_numba"]
            self.use_parallel = config["use_parallel"]
            
            print(f"\nTesting configuration: {config['name']}")
            print(f"  - Numba: {'enabled' if self.use_numba else 'disabled'}")
            print(f"  - Parallel: {'enabled' if self.use_parallel else 'disabled'}")
            
            # Run benchmark
            times = []
            for i in range(n_repeats):
                start_time = time.time()
                traces = self.trace_multiple_lines(start_points_re, time_hours, **kwargs)
                execution_time = time.time() - start_time
                times.append(execution_time)
                print(f"  Run {i+1}/{n_repeats}: {execution_time:.4f} seconds")
            
            # Store results
            avg_time = sum(times) / len(times)
            results[config["name"]] = {
                "times": times,
                "average": avg_time,
                "min": min(times),
                "max": max(times),
                "throughput": len(start_points_re) / avg_time,
                "successful_traces": len(traces) if 'traces' in locals() else 0
            }
            
            print(f"  Average: {avg_time:.4f} seconds "
                  f"({results[config['name']]['throughput']:.2f} lines/sec)")
        
        # Restore original settings
        self.use_numba = orig_numba
        self.use_parallel = orig_parallel
        
        # Print summary
        print("\nBenchmark Summary:")
        baseline = results["baseline"]["average"]
        for name, data in results.items():
            speedup = baseline / data["average"] if name != "baseline" else 1.0
            print(f"  {name}: {data['average']:.4f} sec, "
                  f"throughput: {data['throughput']:.2f} lines/sec, "
                  f"speedup: {speedup:.2f}x")
        
        return results


def plot_vector_field_traces(traces, tracer, title_suffix="", backend='matplotlib',
                            show_earth=True, interactive=True):
    """
    Plot vector field traces with field-type specific styling.
    
    Parameters:
    -----------
    traces : list
        List of trace result dictionaries
    tracer : KamodoVectorFieldTracer
        The tracer object used to generate the traces
    title_suffix : str
        Additional text to add to plot titles
    backend : str
        'matplotlib' or 'plotly' for plotting backend
    show_earth : bool
        Whether to show Earth sphere/circle
    interactive : bool
        For plotly: whether to show interactive features

    Returns:
    --------
    fig : matplotlib.figure.Figure or plotly.graph_objects.Figure
        The generated figure
    """
    if backend.lower() == 'plotly':
        return _plot_vector_field_traces_plotly(
            traces, tracer, title_suffix, show_earth, interactive)
    else:
        return _plot_vector_field_traces_matplotlib(
            traces, tracer, title_suffix, show_earth)


def _plot_vector_field_traces_matplotlib(traces, tracer, title_suffix="", show_earth=True):
    """
    Create matplotlib version of vector field trace plots with aspect ratio 1.
    
    Parameters:
    -----------
    traces : list
        List of trace result dictionaries
    tracer : KamodoVectorFieldTracer
        The tracer object used to generate the traces
    title_suffix : str
        Additional text to add to plot titles
    show_earth : bool
        Whether to show Earth sphere/circle
        
    Returns:
    --------
    fig : matplotlib.figure.Figure
        The generated figure
    """
    fig = plt.figure(figsize=(16, 12))

    # Create subplots
    ax1 = fig.add_subplot(221, projection='3d')  # 3D view
    ax2 = fig.add_subplot(222)  # XY plane
    ax3 = fig.add_subplot(223)  # XZ plane
    ax4 = fig.add_subplot(224)  # Field strength along trace
    
    # Set up colormap based on field type
    if tracer.field_type == 'magnetic':
        cmap = get_colormap('viridis')
        field_unit = 'nT'
    elif tracer.field_type == 'velocity':
        cmap = get_colormap('plasma')
        field_unit = 'km/s'
    elif tracer.field_type == 'current':
        cmap = get_colormap('inferno')
        field_unit = 'A/m²'
    elif tracer.field_type == 'electric':
        cmap = get_colormap('magma')
        field_unit = 'mV/m'
    else:
        cmap = get_colormap('viridis')
        field_unit = 'unit'
    
    # Draw Earth
    if show_earth:
        r_earth = 1.0  # Earth radius in RE
        
        # Earth in 3D view
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = r_earth * np.outer(np.cos(u), np.sin(v))
        y = r_earth * np.outer(np.sin(u), np.sin(v))
        z = r_earth * np.outer(np.ones(np.size(u)), np.cos(v))
        ax1.plot_surface(x, y, z, color='lightblue', alpha=0.3)
        
        # Earth in XY plane
        circle = plt.Circle((0, 0), r_earth, color='lightblue', alpha=0.3)
        ax2.add_patch(circle)
        
        # Earth in XZ plane
        circle = plt.Circle((0, 0), r_earth, color='lightblue', alpha=0.3)
        ax3.add_patch(circle)
    
    # Collect all field magnitudes for consistent color scaling
    all_field_mags = []
    for trace in traces:
        if 'combined' in trace:
            all_field_mags.extend(trace['combined']['field_magnitude'])
        elif 'forward' in trace:
            all_field_mags.extend(trace['forward']['field_magnitude'])
    
    if all_field_mags:
        vmin, vmax = np.min(all_field_mags), np.max(all_field_mags)
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = plt.Normalize(vmin=0, vmax=1)
    
    # Plot each trace
    for i, trace in enumerate(traces):
        if 'combined' in trace:
            coords = trace['combined']['coordinates']
            field_mags = trace['combined']['field_magnitude']
        elif 'forward' in trace:
            coords = trace['forward']['coordinates']
            field_mags = trace['forward']['field_magnitude']
        else:
            continue
        
        # Colors based on field magnitude
        colors = cmap(norm(field_mags))
        
        # 3D plot
        ax1.plot(coords[:, 0], coords[:, 1], coords[:, 2], lw=1.5, alpha=0.7)
        
        # XY plane - plot segments with colors
        for j in range(len(coords)-1):
            ax2.plot(coords[j:j+2, 0], coords[j:j+2, 1], c=colors[j], lw=1.5, alpha=0.7)
        
        # XZ plane - plot segments with colors
        for j in range(len(coords)-1):
            ax3.plot(coords[j:j+2, 0], coords[j:j+2, 2], c=colors[j], lw=1.5, alpha=0.7)
        
        # Field magnitude along trace
        arc_length = np.zeros(len(coords))
        for j in range(1, len(coords)):
            arc_length[j] = arc_length[j-1] + np.linalg.norm(coords[j] - coords[j-1])
        
        ax4.plot(arc_length, field_mags, lw=1.5, label=f"Line {i+1}" if len(traces) <= 10 else "")
    
    # Set titles and labels
    field_type_label = tracer.field_type.capitalize()
    ax1.set_title(f"3D {field_type_label} Field Lines {title_suffix}")
    ax2.set_title(f"XY Plane {title_suffix}")
    ax3.set_title(f"XZ Plane {title_suffix}")
    ax4.set_title(f"{field_type_label} Field Magnitude Along Field Lines")
    
    ax1.set_xlabel('X (RE)')
    ax1.set_ylabel('Y (RE)')
    ax1.set_zlabel('Z (RE)')
    
    ax2.set_xlabel('X (RE)')
    ax2.set_ylabel('Y (RE)')
    
    ax3.set_xlabel('X (RE)')
    ax3.set_ylabel('Z (RE)')
    
    ax4.set_xlabel('Arc Length (RE)')
    ax4.set_ylabel(f'Field Magnitude ({field_unit})')
    ax4.grid(True)
    
    # Equal aspect ratio for spatial plots
    try:
        ax1.set_box_aspect([1,1,1])
    except AttributeError:
        pass  # Older matplotlib versions don't have this method
    ax2.set_aspect('equal')
    ax3.set_aspect('equal')
    
    # Add legend to field magnitude plot
    if len(traces) <= 10:  # Only show legend for reasonable number of lines
        ax4.legend()
    
    # Add colorbar
    if all_field_mags:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=[ax1, ax2, ax3], orientation='vertical', pad=0.01)
        cbar.set_label(f'{field_type_label} Field Magnitude ({field_unit})')
    
    # Adjust layout
    plt.tight_layout()
    
    return fig


def _plot_vector_field_traces_plotly(traces, tracer, title_suffix="", show_earth=True, interactive=True):
    """
    Create plotly version of vector field trace plots.
    Fixed to handle numpy arrays properly in plotly.
    
    Parameters:
    -----------
    traces : list
        List of trace result dictionaries
    tracer : KamodoVectorFieldTracer
        The tracer object used to generate the traces
    title_suffix : str
        Additional text to add to plot titles
    show_earth : bool
        Whether to show Earth sphere/circle
    interactive : bool
        Whether to include interactive elements
        
    Returns:
    --------
    fig : plotly.graph_objects.Figure
        The generated figure
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        import plotly.colors as pc
    except ImportError:
        print("Plotly is not installed. Please install with: pip install plotly")
        return None
    
    # Create figure with 2x2 subplots
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{'type': 'scene'}, {'type': 'xy'}],
               [{'type': 'xy'}, {'type': 'xy'}]],
        subplot_titles=(
            f"3D {tracer.field_type.capitalize()} Field Lines {title_suffix}",
            f"XY Plane {title_suffix}",
            f"XZ Plane {title_suffix}",
            f"{tracer.field_type.capitalize()} Field Magnitude Along Field Lines"
        )
    )
    
    # Set up color based on field type
    if tracer.field_type == 'magnetic':
        colorscale = 'Viridis'
        field_unit = 'nT'
    elif tracer.field_type == 'velocity':
        colorscale = 'Plasma'
        field_unit = 'km/s'
    elif tracer.field_type == 'current':
        colorscale = 'Inferno'
        field_unit = 'A/m²'
    elif tracer.field_type == 'electric':
        colorscale = 'Magma'
        field_unit = 'mV/m'
    else:
        colorscale = 'Viridis'
        field_unit = 'unit'
    
    # Draw Earth
    if show_earth:
        r_earth = 1.0  # Earth radius in RE
        
        # Create a sphere for Earth in 3D view
        theta = np.linspace(0, 2 * np.pi, 50)
        phi = np.linspace(0, np.pi, 50)
        x = r_earth * np.outer(np.cos(theta), np.sin(phi))
        y = r_earth * np.outer(np.sin(theta), np.sin(phi))
        z = r_earth * np.outer(np.ones(np.size(theta)), np.cos(phi))
        
        fig.add_trace(
            go.Surface(
                x=x, y=y, z=z,
                colorscale=[[0, 'lightblue'], [1, 'lightblue']],
                opacity=0.3,
                showscale=False,
                name='Earth'
            ),
            row=1, col=1
        )
        
        # Add circles for Earth in 2D views
        theta_circle = np.linspace(0, 2*np.pi, 100)
        x_circle = r_earth * np.cos(theta_circle)
        y_circle = r_earth * np.sin(theta_circle)
        
        # XY plane (Earth)
        fig.add_trace(
            go.Scatter(
                x=x_circle, y=y_circle,
                mode='lines',
                line=dict(color='lightblue', width=2),
                fill='toself',
                fillcolor='rgba(173, 216, 230, 0.3)',
                name='Earth',
                showlegend=False
            ),
            row=1, col=2
        )
        
        # XZ plane (Earth)
        fig.add_trace(
            go.Scatter(
                x=x_circle, y=np.zeros(100),
                mode='lines',
                line=dict(color='lightblue', width=2),
                fill='toself',
                fillcolor='rgba(173, 216, 230, 0.3)',
                showlegend=False
            ),
            row=2, col=1
        )
    
    # Calculate global field magnitude range for consistent coloring
    all_field_mags = []
    for trace in traces:
        if 'combined' in trace:
            all_field_mags.extend(trace['combined']['field_magnitude'])
        elif 'forward' in trace:
            all_field_mags.extend(trace['forward']['field_magnitude'])
    
    if all_field_mags:
        field_min, field_max = np.min(all_field_mags), np.max(all_field_mags)
    else:
        field_min, field_max = 0, 1
    
    # Plot each trace
    for i, trace in enumerate(traces):
        if 'combined' in trace:
            coords = trace['combined']['coordinates']
            field_mags = trace['combined']['field_magnitude']
        elif 'forward' in trace:
            coords = trace['forward']['coordinates']
            field_mags = trace['forward']['field_magnitude']
        else:
            continue
        
        # 3D plot - this works with arrays for color
        fig.add_trace(
            go.Scatter3d(
                x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
                mode='lines',
                line=dict(
                    color=field_mags,
                    colorscale=colorscale,
                    width=4,
                    cmin=field_min,
                    cmax=field_max,
                    colorbar=dict(
                        title=f'{field_unit}',
                        thickness=20,
                        len=0.7,
                        x=1.02
                    ) if i == 0 else None,  # Only show colorbar for first trace
                ),
                name=f'Line {i+1}',
                showlegend=False
            ),
            row=1, col=1
        )
        
        # For 2D plots, use average color approach to avoid numpy array issues
        avg_field_mag = np.mean(field_mags)
        
        # Normalize to 0-1 range
        if field_max > field_min:
            normalized_mag = (avg_field_mag - field_min) / (field_max - field_min)
        else:
            normalized_mag = 0.5
        
        # Get color from colorscale
        colors = pc.sample_colorscale(colorscale, [normalized_mag])
        line_color = colors[0]
        
        # XY plane
        fig.add_trace(
            go.Scatter(
                x=coords[:, 0], y=coords[:, 1],
                mode='lines',
                line=dict(
                    color=line_color,
                    width=3,
                ),
                name=f'Line {i+1}' if len(traces) <= 10 else "",
                showlegend=False,
                hovertemplate=f'X: %{{x:.2f}} RE<br>Y: %{{y:.2f}} RE<br>Avg Field: {avg_field_mag:.2f} {field_unit}<extra></extra>'
            ),
            row=1, col=2
        )
        
        # XZ plane
        fig.add_trace(
            go.Scatter(
                x=coords[:, 0], y=coords[:, 2],
                mode='lines',
                line=dict(
                    color=line_color,
                    width=3,
                ),
                showlegend=False,
                hovertemplate=f'X: %{{x:.2f}} RE<br>Z: %{{y:.2f}} RE<br>Avg Field: {avg_field_mag:.2f} {field_unit}<extra></extra>'
            ),
            row=2, col=1
        )
        
        # Field magnitude along trace
        arc_length = np.zeros(len(coords))
        for j in range(1, len(coords)):
            arc_length[j] = arc_length[j-1] + np.linalg.norm(coords[j] - coords[j-1])
        
        fig.add_trace(
            go.Scatter(
                x=arc_length, y=field_mags,
                mode='lines',
                line=dict(color=line_color, width=2),
                name=f'Line {i+1}' if len(traces) <= 10 else "",
                showlegend=len(traces) <= 10,
                hovertemplate='Arc Length: %{x:.2f} RE<br>Field Magnitude: %{y:.2f} ' + field_unit + '<extra></extra>'
            ),
            row=2, col=2
        )
    
    # Update axis labels
    fig.update_scenes(
        xaxis_title='X (RE)',
        yaxis_title='Y (RE)',
        zaxis_title='Z (RE)',
        aspectmode='data'
    )
    
    fig.update_xaxes(title_text='X (RE)', row=1, col=2)
    fig.update_yaxes(title_text='Y (RE)', row=1, col=2)
    
    fig.update_xaxes(title_text='X (RE)', row=2, col=1)
    fig.update_yaxes(title_text='Z (RE)', row=2, col=1)
    
    fig.update_xaxes(title_text='Arc Length (RE)', row=2, col=2)
    fig.update_yaxes(title_text=f'Field Magnitude ({field_unit})', row=2, col=2)
    
    # Equal aspect ratio for 2D plots
    fig.update_xaxes(scaleanchor='y', scaleratio=1, row=1, col=2)
    fig.update_xaxes(scaleanchor='y', scaleratio=1, row=2, col=1)
    
    # Update layout
    fig.update_layout(
        title=f'{tracer.field_type.capitalize()} Field Line Tracing {title_suffix}',
        width=1200,
        height=900,
        hovermode='closest' if interactive else False
    )
    
    return fig


def create_magnetic_tracer(kamodo_object, component_names=None, use_numba=True, use_parallel=False, n_jobs=-1, verbose=False):
    """
    Create a tracer specifically for magnetic field lines.
    
    Parameters:
    -----------
    kamodo_object : object
        Kamodo object with magnetic field components
    component_names : list, optional
        List of field component variable names in kamodo object.
        Default: ['B_x', 'B_y', 'B_z']
    use_numba : bool, optional (default=True)
        Whether to use numba optimization
    use_parallel : bool, optional (default=False)
        Whether to use parallel processing
    n_jobs : int, optional (default=-1)
        Number of jobs for parallel processing
    verbose : bool, optional (default=False)
        Whether to print detailed information
    
    Returns:
    --------
    tracer : KamodoVectorFieldTracer
        Configured tracer for magnetic field lines
    """
    if component_names is None:
        component_names = ['B_x', 'B_y', 'B_z']
    
    return KamodoVectorFieldTracer(
        kamodo_object, 
        field_type='magnetic', 
        component_names=component_names,
        use_numba=use_numba,
        use_parallel=use_parallel,
        n_jobs=n_jobs,
        verbose=verbose
    )


def create_velocity_tracer(kamodo_object, component_names=None, use_numba=True, use_parallel=False, n_jobs=-1, verbose=False):
    """
    Create a tracer specifically for velocity streamlines.
    
    Parameters:
    -----------
    kamodo_object : object
        Kamodo object with velocity field components
    component_names : list, optional
        List of field component variable names in kamodo object.
        Default: ['v_x', 'v_y', 'v_z']
    use_numba : bool, optional (default=True)
        Whether to use numba optimization
    use_parallel : bool, optional (default=False)
        Whether to use parallel processing
    n_jobs : int, optional (default=-1)
        Number of jobs for parallel processing
    verbose : bool, optional (default=False)
        Whether to print detailed information
    
    Returns:
    --------
    tracer : KamodoVectorFieldTracer
        Configured tracer for velocity streamlines
    """
    if component_names is None:
        component_names = ['v_x', 'v_y', 'v_z']
    
    return KamodoVectorFieldTracer(
        kamodo_object, 
        field_type='velocity', 
        component_names=component_names,
        use_numba=use_numba,
        use_parallel=use_parallel,
        n_jobs=n_jobs,
        verbose=verbose
    )


def create_current_tracer(kamodo_object, component_names=None, use_numba=True, use_parallel=False, n_jobs=-1, verbose=False):
    """
    Create a tracer specifically for current density streamlines.
    
    Parameters:
    -----------
    kamodo_object : object
        Kamodo object with current density field components
    component_names : list, optional
        List of field component variable names in kamodo object.
        Default: ['J_x', 'J_y', 'J_z']
    use_numba : bool, optional (default=True)
        Whether to use numba optimization
    use_parallel : bool, optional (default=False)
        Whether to use parallel processing
    n_jobs : int, optional (default=-1)
        Number of jobs for parallel processing
    verbose : bool, optional (default=False)
        Whether to print detailed information
    
    Returns:
    --------
    tracer : KamodoVectorFieldTracer
        Configured tracer for current streamlines
    """
    if component_names is None:
        component_names = ['J_x', 'J_y', 'J_z']
    
    return KamodoVectorFieldTracer(
        kamodo_object, 
        field_type='current', 
        component_names=component_names,
        use_numba=use_numba,
        use_parallel=use_parallel,
        n_jobs=n_jobs,
        verbose=verbose
    )


def create_electric_tracer(kamodo_object, component_names=None, use_numba=True, use_parallel=False, n_jobs=-1, verbose=False):
    """
    Create a tracer specifically for electric field lines.
    
    Parameters:
    -----------
    kamodo_object : object
        Kamodo object with electric field components
    component_names : list, optional
        List of field component variable names in kamodo object.
        Default: ['E_x', 'E_y', 'E_z']
    use_numba : bool, optional (default=True)
        Whether to use numba optimization
    use_parallel : bool, optional (default=False)
        Whether to use parallel processing
    n_jobs : int, optional (default=-1)
        Number of jobs for parallel processing
    verbose : bool, optional (default=False)
        Whether to print detailed information
    
    Returns:
    --------
    tracer : KamodoVectorFieldTracer
        Configured tracer for electric field lines
    """
    if component_names is None:
        component_names = ['E_x', 'E_y', 'E_z']
    
    return KamodoVectorFieldTracer(
        kamodo_object, 
        field_type='electric', 
        component_names=component_names,
        use_numba=use_numba,
        use_parallel=use_parallel,
        n_jobs=n_jobs,
        verbose=verbose
    )


def inspect_kamodo_object(kamodo_object, verbose=True):
    """
    Inspect Kamodo object for vector field components and their calling signatures.
    
    Parameters:
    -----------
    kamodo_object : object
        Kamodo object to inspect
    verbose : bool
        Whether to print detailed information
    
    Returns:
    --------
    field_info : dict
        Dictionary of detected field types and their detailed information
    """
    if not hasattr(kamodo_object, 'keys'):
        if verbose:
            print("Error: Object does not appear to be a valid Kamodo object")
        return {}
        
    variables = list(kamodo_object.keys())
    if verbose:
        print(f"Found {len(variables)} variables in Kamodo object:")
    
    # Inspect function signatures
    function_info = {}
    for var_name in variables:
        var_func = kamodo_object[var_name]
        if hasattr(var_func, '__call__'):
            try:
                import inspect
                sig = inspect.signature(var_func)
                param_names = list(sig.parameters.keys())
                function_info[var_name] = {
                    'parameters': param_names,
                    'signature': str(sig)
                }
                if verbose:
                    print(f"  {var_name}: {sig}")
            except Exception as e:
                if verbose:
                    print(f"  {var_name}: Could not inspect signature - {e}")
                function_info[var_name] = {'error': str(e)}
        else:
            if verbose:
                print(f"  {var_name}: Not callable (constant?)")
            function_info[var_name] = {'type': 'constant'}
    
    # Common field component patterns
    field_patterns = {
        'magnetic': [['B_x', 'B_y', 'B_z'], ['b_x', 'b_y', 'b_z'], ['Bx', 'By', 'Bz']],
        'velocity': [['v_x', 'v_y', 'v_z'], ['u_x', 'u_y', 'u_z'], ['vx', 'vy', 'vz']],
        'electric': [['E_x', 'E_y', 'E_z'], ['e_x', 'e_y', 'e_z'], ['Ex', 'Ey', 'Ez']],
        'current': [['J_x', 'J_y', 'J_z'], ['j_x', 'j_y', 'j_z'], ['Jx', 'Jy', 'Jz']],
    }
    
    # Check which fields are available
    available_fields = {}
    for field_type, pattern_lists in field_patterns.items():
        for patterns in pattern_lists:
            components = [var for var in variables if var in patterns]
            if len(components) == 3:  # Must have all three components
                if verbose:
                    print(f"✓ {field_type.capitalize()} field components found: {components}")
                
                # Get signature info for these components
                component_sigs = {}
                for comp in components:
                    if comp in function_info:
                        component_sigs[comp] = function_info[comp]
                
                available_fields[field_type] = {
                    'components': components,
                    'signatures': component_sigs
                }
                break
        else:
            if verbose:
                print(f"✗ No complete {field_type} field components found")
    
    if verbose:
        print(f"\nAll available variables: {variables}")
        print("\nKamodo calling format: component([time, x, y, z])")
    
    return {
        'available_fields': available_fields,
        'all_functions': function_info,
        'all_variables': variables
    }


def benchmark_vectorfield_tracer(kamodo_object, n_lines=10, use_numba_options=[True, False],
                               use_parallel_options=[True, False], field_type='magnetic'):
    """
    Comprehensive benchmark of vector field tracer with different optimization settings.
    
    Parameters:
    -----------
    kamodo_object : object
        Kamodo object with field components
    n_lines : int
        Number of field lines to trace for benchmark
    use_numba_options : list
        List of boolean options to test for numba usage
    use_parallel_options : list
        List of boolean options to test for parallel processing
    field_type : str
        Type of field to benchmark ('magnetic', 'velocity', etc.)
        
    Returns:
    --------
    results : dict
        Dictionary of benchmark results
    """
    print(f"Starting comprehensive benchmark with {n_lines} field lines")
    
    # Generate random start points in magnetosphere
    np.random.seed(42)  # For reproducibility
    theta = np.random.uniform(0, 2*np.pi, n_lines)
    phi = np.random.uniform(0, np.pi, n_lines)
    r = np.random.uniform(3, 10, n_lines)
    
    start_points = np.array([
        r * np.cos(theta) * np.sin(phi),
        r * np.sin(theta) * np.sin(phi),
        r * np.cos(phi)
    ]).T
    
    results = {}
    
    # Detect available fields in kamodo object
    field_info = inspect_kamodo_object(kamodo_object, verbose=False)
    available_fields = field_info['available_fields']
    
    if field_type not in available_fields:
        print(f"Error: {field_type} field not found in Kamodo object")
        return {}
    
    components = available_fields[field_type]['components']
    print(f"\nBenchmarking {field_type} field tracing with components: {components}")
    field_results = {}
    
    for use_numba in use_numba_options:
        for use_parallel in use_parallel_options:
            # Skip combinations that aren't available
            if use_numba and not NUMBA_AVAILABLE:
                continue
            if use_parallel and not JOBLIB_AVAILABLE:
                continue
                
            config_name = f"numba={use_numba}_parallel={use_parallel}"
            print(f"\nConfiguration: {config_name}")
            
            # Create tracer with specified settings
            tracer = KamodoVectorFieldTracer(
                kamodo_object,
                field_type=field_type,
                component_names=components,
                use_numba=use_numba,
                use_parallel=use_parallel,
                verbose=False
            )
            
            # Run benchmark
            start_time = time.time()
            traces = tracer.trace_multiple_lines(
                start_points, time_hours=0.0, 
                max_steps=1000, step_size_re=0.1
            )
            end_time = time.time()
            
            # Calculate metrics
            execution_time = end_time - start_time
            lines_per_second = n_lines / execution_time if execution_time > 0 else 0
            avg_points_per_line = np.mean([len(t['combined']['coordinates']) 
                                           if 'combined' in t else 
                                           len(t['forward']['coordinates']) 
                                           for t in traces]) if traces else 0
            
            # Store results
            field_results[config_name] = {
                'execution_time': execution_time,
                'lines_per_second': lines_per_second,
                'avg_points_per_line': avg_points_per_line,
                'num_traces': len(traces)
            }
            
            print(f"Execution time: {execution_time:.4f} seconds")
            print(f"Performance: {lines_per_second:.2f} lines/second")
            print(f"Avg points per line: {avg_points_per_line:.1f}")
            print(f"Successful traces: {len(traces)}/{n_lines}")
    
    # Store results for this field type
    results[field_type] = field_results
    
    # Print overall summary
    print(f"\n===== BENCHMARK SUMMARY FOR {field_type.upper()} FIELD =====")
    
    # Find baseline (no optimization) config
    baseline = field_results.get('numba=False_parallel=False', None)
    if baseline is None:
        print("  Baseline configuration not found")
        return results
    
    baseline_time = baseline['execution_time']
    
    # Compare each configuration
    for config_name, config_results in field_results.items():
        speedup = baseline_time / config_results['execution_time'] if config_results['execution_time'] > 0 else 0
        print(f"  {config_name}: {config_results['execution_time']:.4f} sec, "
              f"{speedup:.2f}x speedup, "
              f"{config_results['lines_per_second']:.2f} lines/sec")
    
    return results


def memory_efficient_trace(tracer, start_points, time_hours, max_memory_mb=1000, **kwargs):
    """
    Memory-efficient field line tracing for very large datasets.
    
    This function automatically determines batch sizes based on available memory
    and processes start points in chunks to avoid memory overflow.
    
    Parameters:
    -----------
    tracer : KamodoVectorFieldTracer
        The tracer object to use
    start_points : array_like
        Array of shape (n_points, 3) containing start positions
    time_hours : float
        Time in hours
    max_memory_mb : float
        Maximum memory usage in MB (approximate)
    **kwargs : dict
        Additional arguments passed to trace_multiple_lines
        
    Returns:
    --------
    results : list
        List of all trace results
    memory_stats : dict
        Dictionary with memory usage statistics
    """
    try:
        import psutil
        process = psutil.Process()
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        peak_memory = initial_memory
        memory_available = True
    except ImportError:
        print("psutil not available. Memory monitoring disabled.")
        memory_available = False
        initial_memory = 0
        peak_memory = 0
    
    start_points = np.array(start_points)
    n_points = len(start_points)
    
    # Estimate memory per trace (rough approximation)
    # Assume ~1000 points per trace, 3 coords + 1 magnitude + 3 components = 7 floats
    # Plus overhead, estimate ~10KB per trace
    memory_per_trace_kb = 10
    max_traces_per_batch = int((max_memory_mb * 1000) / memory_per_trace_kb)
    batch_size = min(max_traces_per_batch, 1000)  # Cap at 1000 for reasonable processing
    
    print(f"Processing {n_points} points with estimated batch size: {batch_size}")
    print(f"Target memory limit: {max_memory_mb} MB")
    
    all_results = []
    n_batches = (n_points + batch_size - 1) // batch_size
    
    for i in range(n_batches):
        start_idx = i * batch_size
        end_idx = min((i + 1) * batch_size, n_points)
        batch = start_points[start_idx:end_idx]
        
        print(f"Processing batch {i+1}/{n_batches}: {len(batch)} points")
        
        # Process batch
        batch_results = tracer.trace_multiple_lines(batch, time_hours, **kwargs)
        all_results.extend(batch_results)
        
        # Check memory usage
        if memory_available:
            current_memory = process.memory_info().rss / 1024 / 1024
            peak_memory = max(peak_memory, current_memory)
            print(f"  Memory usage: {current_memory:.1f} MB (peak: {peak_memory:.1f} MB)")
        
        # Force garbage collection every few batches
        if i % 3 == 0:
            gc.collect()
            if tracer.verbose:
                print("  Performed garbage collection")
    
    # Final memory stats
    final_memory = process.memory_info().rss / 1024 / 1024 if memory_available else 0
    memory_stats = {
        'initial_memory_mb': initial_memory,
        'peak_memory_mb': peak_memory,
        'final_memory_mb': final_memory,
        'memory_increase_mb': final_memory - initial_memory,
        'batches_processed': n_batches,
        'points_per_batch': batch_size,
        'memory_monitoring_available': memory_available
    }
    
    if memory_available:
        print(f"\nMemory Summary:")
        print(f"  Initial: {initial_memory:.1f} MB")
        print(f"  Peak: {peak_memory:.1f} MB")
        print(f"  Final: {final_memory:.1f} MB")
        print(f"  Increase: {memory_stats['memory_increase_mb']:.1f} MB")
    
    return all_results, memory_stats


def save_traces_to_file(traces, filename, format='npz', include_metadata=True):
    """
    Save field line traces to file for later analysis.
    
    Parameters:
    -----------
    traces : list
        List of trace result dictionaries
    filename : str
        Output filename
    format : str
        Output format ('npz', 'hdf5', 'json')
    include_metadata : bool
        Whether to include metadata like execution times
        
    Returns:
    --------
    success : bool
        Whether save was successful
    """
    try:
        if format.lower() == 'npz':
            # Save as compressed numpy arrays
            save_dict = {}
            for i, trace in enumerate(traces):
                prefix = f'trace_{i:04d}_'
                
                if 'combined' in trace:
                    save_dict[f'{prefix}coords'] = trace['combined']['coordinates']
                    save_dict[f'{prefix}field_mag'] = trace['combined']['field_magnitude']
                    save_dict[f'{prefix}field_comp'] = trace['combined']['field_components']
                    save_dict[f'{prefix}length'] = trace['combined']['length']
                elif 'forward' in trace:
                    save_dict[f'{prefix}coords'] = trace['forward']['coordinates']
                    save_dict[f'{prefix}field_mag'] = trace['forward']['field_magnitude']
                    save_dict[f'{prefix}field_comp'] = trace['forward']['field_components']
                    save_dict[f'{prefix}length'] = trace['forward']['length']
                
                if include_metadata and 'execution_time' in trace:
                    save_dict[f'{prefix}exec_time'] = trace['execution_time']
            
            np.savez_compressed(filename, **save_dict)
            print(f"Saved {len(traces)} traces to {filename}")
            
        elif format.lower() == 'hdf5':
            try:
                import h5py
                with h5py.File(filename, 'w') as f:
                    for i, trace in enumerate(traces):
                        grp = f.create_group(f'trace_{i:04d}')
                        
                        if 'combined' in trace:
                            data = trace['combined']
                        elif 'forward' in trace:
                            data = trace['forward']
                        else:
                            continue
                        
                        grp.create_dataset('coordinates', data=data['coordinates'])
                        grp.create_dataset('field_magnitude', data=data['field_magnitude'])
                        grp.create_dataset('field_components', data=data['field_components'])
                        grp.create_dataset('length', data=data['length'])
                        
                        if include_metadata and 'execution_time' in trace:
                            grp.attrs['execution_time'] = trace['execution_time']
                
                print(f"Saved {len(traces)} traces to {filename} (HDF5)")
                
            except ImportError:
                print("h5py not available. Install with: pip install h5py")
                return False
                
        elif format.lower() == 'json':
            import json
            
            # Convert numpy arrays to lists for JSON serialization
            json_traces = []
            for trace in traces:
                json_trace = {}
                
                if 'combined' in trace:
                    data = trace['combined']
                elif 'forward' in trace:
                    data = trace['forward']
                else:
                    continue
                
                json_trace['coordinates'] = data['coordinates'].tolist()
                json_trace['field_magnitude'] = data['field_magnitude'].tolist()
                json_trace['field_components'] = data['field_components'].tolist()
                json_trace['length'] = float(data['length'])
                json_trace['n_points'] = int(data['n_points'])
                
                if include_metadata and 'execution_time' in trace:
                    json_trace['execution_time'] = trace['execution_time']
                
                json_traces.append(json_trace)
            
            with open(filename, 'w') as f:
                json.dump(json_traces, f, indent=2)
            
            print(f"Saved {len(traces)} traces to {filename} (JSON)")
        
        else:
            print(f"Unsupported format: {format}")
            return False
        
        return True
        
    except Exception as e:
        print(f"Error saving traces: {e}")
        return False


def load_traces_from_file(filename, format='npz'):
    """
    Load field line traces from file.
    
    Parameters:
    -----------
    filename : str
        Input filename
    format : str
        Input format ('npz', 'hdf5', 'json')
        
    Returns:
    --------
    traces : list
        List of loaded trace dictionaries
    """
    try:
        if format.lower() == 'npz':
            data = np.load(filename)
            traces = []
            
            # Group data by trace number
            trace_nums = set()
            for key in data.keys():
                if '_coords' in key:
                    trace_num = key.split('_coords')[0]
                    trace_nums.add(trace_num)
            
            for trace_num in sorted(trace_nums):
                trace = {
                    'combined': {
                        'coordinates': data[f'{trace_num}_coords'],
                        'field_magnitude': data[f'{trace_num}_field_mag'],
                        'field_components': data[f'{trace_num}_field_comp'],
                        'length': float(data[f'{trace_num}_length']),
                        'n_points': len(data[f'{trace_num}_coords'])
                    }
                }
                
                if f'{trace_num}_exec_time' in data:
                    trace['execution_time'] = float(data[f'{trace_num}_exec_time'])
                
                traces.append(trace)
            
            print(f"Loaded {len(traces)} traces from {filename}")
            
        elif format.lower() == 'hdf5':
            try:
                import h5py
                traces = []
                
                with h5py.File(filename, 'r') as f:
                    for trace_name in sorted(f.keys()):
                        grp = f[trace_name]
                        
                        trace = {
                            'combined': {
                                'coordinates': np.array(grp['coordinates']),
                                'field_magnitude': np.array(grp['field_magnitude']),
                                'field_components': np.array(grp['field_components']),
                                'length': float(grp['length'][()]),
                                'n_points': len(grp['coordinates'])
                            }
                        }
                        
                        if 'execution_time' in grp.attrs:
                            trace['execution_time'] = grp.attrs['execution_time']
                        
                        traces.append(trace)
                
                print(f"Loaded {len(traces)} traces from {filename} (HDF5)")
                
            except ImportError:
                print("h5py not available. Install with: pip install h5py")
                return []
                
        elif format.lower() == 'json':
            import json
            
            with open(filename, 'r') as f:
                json_traces = json.load(f)
            
            traces = []
            for json_trace in json_traces:
                trace = {
                    'combined': {
                        'coordinates': np.array(json_trace['coordinates']),
                        'field_magnitude': np.array(json_trace['field_magnitude']),
                        'field_components': np.array(json_trace['field_components']),
                        'length': json_trace['length'],
                        'n_points': json_trace['n_points']
                    }
                }
                
                if 'execution_time' in json_trace:
                    trace['execution_time'] = json_trace['execution_time']
                
                traces.append(trace)
            
            print(f"Loaded {len(traces)} traces from {filename} (JSON)")
        
        else:
            print(f"Unsupported format: {format}")
            return []
        
        return traces
        
    except Exception as e:
        print(f"Error loading traces: {e}")
        return []


def test_kamodo_integration(kamodo_object, component_names=['B_x', 'B_y', 'B_z'], 
                          test_time=2.0, test_position=[10, 0, 0]):
    """
    Test Kamodo integration with the correct calling convention.
    
    Parameters:
    -----------
    kamodo_object : object
        Kamodo object to test
    component_names : list
        List of component names to test
    test_time : float
        Test time in hours
    test_position : list
        Test position [x, y, z] in RE
    """
    print("Testing Kamodo Integration")
    print("=" * 30)
    
    test_coords = [test_time] + test_position
    print(f"Test coordinates: [t={test_time}, x={test_position[0]}, y={test_position[1]}, z={test_position[2]}]")
    
    for comp_name in component_names:
        if comp_name in kamodo_object:
            try:
                # Test single point call
                result = kamodo_object[comp_name](test_coords)
                print(f"✓ {comp_name}({test_coords}) = {result}")
                
                # Test multiple point call
                test_coords2 = [test_time, test_position[0]+1, test_position[1], test_position[2]]
                multi_result = kamodo_object[comp_name]([test_coords, test_coords2])
                print(f"✓ {comp_name} multi-point: {multi_result}")
                
            except Exception as e:
                print(f"✗ Error calling {comp_name}: {e}")
        else:
            print(f"✗ Component {comp_name} not found in Kamodo object")
    
    print("\nTest completed.")


def example_usage():
    """
    Example showing how to use the optimized vector field tracer.
    
    This function demonstrates the main features and usage patterns
    of the optimized vector field tracer.
    """
    print("=== Optimized Vector Field Tracer Usage Examples ===\n")
    
    print("1. Basic Usage - Create a magnetic field tracer:")
    print("   tracer = create_magnetic_tracer(kamodo_object)")
    print("   # Automatically detects B_x, B_y, B_z components")
    print("   # Enables numba and parallel processing by default\n")
    
    print("2. Test Kamodo integration:")
    print("   test_kamodo_integration(kamodo_object)")
    print("   # Verifies that Kamodo calling convention works\n")
    
    print("3. Trace a single field line:")
    print("   result = tracer.trace_vector_line(")
    print("       start_position_re=[10, 0, 0],  # Start at [10, 0, 0] RE")
    print("       time_hours=2.0,                # At t=2 hours")
    print("       direction='both',              # Trace both directions")
    print("       adaptive_stepping=True         # Use adaptive step size")
    print("   )")
    print("   # Returns dict with 'forward', 'backward', 'combined' results\n")
    
    print("4. Trace multiple field lines (with automatic batching):")
    print("   start_points = [[10, 0, 0], [10, 5, 0], [10, -5, 0]]")
    print("   results = tracer.trace_multiple_lines(")
    print("       start_points, time_hours=2.0, show_progress=True")
    print("   )")
    print("   # Automatically uses batch processing for large datasets\n")
    
    print("5. Memory-efficient processing for large datasets:")
    print("   large_results, memory_stats = memory_efficient_trace(")
    print("       tracer, large_start_points, time_hours=2.0, max_memory_mb=2000")
    print("   )")
    print("   # Automatically manages memory usage\n")
    
    print("6. Benchmark performance:")
    print("   benchmark_results = tracer.benchmark(")
    print("       start_points, time_hours=2.0, n_repeats=3")
    print("   )")
    print("   # Tests different optimization combinations\n")
    
    print("7. Monitor performance:")
    print("   stats = tracer.get_performance_stats()")
    print("   print(f'Traces per second: {stats[\"traces_per_second\"]:.2f}')\n")
    
    print("8. Visualize field lines:")
    print("   # Matplotlib (default)")
    print("   fig = plot_vector_field_traces(results, tracer)")
    print("   plt.show()")
    print("   ")
    print("   # Plotly (interactive)")
    print("   fig = plot_vector_field_traces(results, tracer, backend='plotly')")
    print("   fig.show()\n")
    
    print("9. Save and load results:")
    print("   save_traces_to_file(results, 'field_lines.npz')")
    print("   loaded_traces = load_traces_from_file('field_lines.npz')\n")
    
    print("10. Inspect Kamodo object:")
    print("    field_info = inspect_kamodo_object(kamodo_object)")
    print("    # Shows available field components and their signatures\n")
    
    print("Performance Tips:")
    print("- Use numba=True for 10-100x speedup on numerical computations")
    print("- Use use_parallel=True for linear scaling with CPU cores")
    print("- Use adaptive_stepping=True for better accuracy with fewer points")
    print("- Use batch processing for very large datasets (>1000 start points)")
    print("- Monitor memory usage with get_performance_stats()")
    print("- Use plotly backend for interactive 3D visualization")
    print("- Set verbose=True for detailed execution information")
    print("\nFor detailed documentation, see the docstrings of each function and class.")


if __name__ == "__main__":
    """
    Main execution section - runs when script is called directly.
    
    This section demonstrates the capabilities of the optimized vector field tracer
    and can be used for testing and validation.
    """
    print("Optimized Kamodo Vector Field Tracer")
    print("====================================")
    
    # Check available optimizations
    print("\nAvailable optimizations:")
    print(f"  Numba JIT compilation: {'✓' if NUMBA_AVAILABLE else '✗'}")
    print(f"  Parallel processing:   {'✓' if JOBLIB_AVAILABLE else '✗'}")
    print(f"  Progress bars (tqdm):  {'✓' if TQDM_AVAILABLE else '✗'}")
    
    if not NUMBA_AVAILABLE:
        print("\n⚠️  For best performance, install numba:")
        print("    pip install numba")
    
    if not JOBLIB_AVAILABLE:
        print("\n⚠️  For parallel processing, install joblib:")
        print("    pip install joblib")
    
    if not TQDM_AVAILABLE:
        print("\n⚠️  For progress bars, install tqdm:")
        print("    pip install tqdm")
    
    # Show example usage
    print("\n" + "="*50)
    example_usage()
    
    print("\n" + "="*50)
    print("For more examples and detailed documentation,")
    print("see the docstrings and comments in this file.")

