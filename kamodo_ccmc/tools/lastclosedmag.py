"""
Kamodo Last Closed Magnetic Field Lines

Functions for finding the last closed magnetic field lines and magnetopause boundary.
All coordinates are in Earth radii (RE).

Updated to use optimized vector field tracer with:
- Batch processing for efficient handling of multiple field lines
- Parallel processing for improved performance
- Memory management for large datasets
- Enhanced progress tracking and benchmarking
"""

import numpy as np
import matplotlib.pyplot as plt
import time
try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

# Import the optimized vector field tracer
import kamodo_ccmc.tools.vectorfieldtracer as vft


def _is_field_line_closed(trace_result):
    """
    Determine if a field line is closed based on where it terminates.

    Parameters:
    -----------
    trace_result : dict
        Result from trace_vector_line

    Returns:
    --------
    is_closed : bool
        True if field line appears to be closed (both ends hit Earth)
    """
    # Use combined trace if available
    if 'combined' in trace_result:
        stop_reasons = trace_result['combined'].get('stop_reasons', [])
    else:
        # Check individual directions
        stop_reasons = []
        if 'forward' in trace_result:
            stop_reasons.append(trace_result['forward']['stop_reason'])
        if 'backward' in trace_result:
            stop_reasons.append(trace_result['backward']['stop_reason'])

    # Count how many ends hit Earth's surface (minimum altitude)
    earth_hits = sum(1 for reason in stop_reasons if 'min_altitude' in reason or 'Earth_surface' in reason)
    
    # A closed field line should have both ends hitting Earth
    return earth_hits >= 2


def _find_last_closed_radius_bisection(tracer, longitude_deg, time_hours, 
                                     r_min=2.0, r_max=15.0, tolerance=0.1,
                                     max_iterations=20, verbose=False):
    """
    Find the last closed field line radius at a given longitude using bisection method.
    
    Parameters:
    -----------
    tracer : vft.KamodoVectorFieldTracer
        Magnetic field tracer
    longitude_deg : float
        Longitude in degrees
    time_hours : float
        Time in hours
    r_min : float
        Minimum radius to search
    r_max : float
        Maximum radius to search
    tolerance : float
        Tolerance for bisection method
    max_iterations : int
        Maximum number of bisection iterations
    verbose : bool
        Whether to print detailed progress
        
    Returns:
    --------
    last_closed_r : float
        Radius of last closed field line (NaN if not found)
    accuracy : float
        Estimated accuracy of the result
    """
    longitude_rad = np.deg2rad(longitude_deg)
    
    # Test function to check if field line at radius r is closed
    def is_closed_at_radius(r):
        start_pos = [r * np.cos(longitude_rad), r * np.sin(longitude_rad), 0.0]
        try:
            trace_result = tracer.trace_vector_line(
                start_pos, time_hours,
                max_steps=2000,
                step_size_re=0.05,
                direction='both'
            )
            return _is_field_line_closed(trace_result)
        except Exception as e:
            if verbose:
                print(f"    Error tracing at r={r:.2f}: {e}")
            return False
    
    # Check if we have any closed field lines in the range
    if not is_closed_at_radius(r_min):
        if verbose:
            print(f"    No closed field lines found at minimum radius {r_min:.2f}")
        return np.nan, np.nan
    
    # Check if all field lines are closed up to maximum radius
    if is_closed_at_radius(r_max):
        if verbose:
            print(f"    All field lines closed up to maximum radius {r_max:.2f}")
        return r_max, tolerance
    
    # Bisection method
    r_low = r_min
    r_high = r_max
    
    for iteration in range(max_iterations):
        r_mid = (r_low + r_high) / 2.0
        
        if verbose:
            print(f"    Iteration {iteration+1}: testing r={r_mid:.3f} RE")
        
        if is_closed_at_radius(r_mid):
            r_low = r_mid  # Last closed is at least r_mid
        else:
            r_high = r_mid  # Last closed is less than r_mid
        
        # Check convergence
        if (r_high - r_low) < tolerance:
            break
    
    last_closed_r = (r_low + r_high) / 2.0
    accuracy = (r_high - r_low) / 2.0
    
    if verbose:
        print(f"    Found last closed at r={last_closed_r:.3f} ± {accuracy:.3f} RE")
    
    return last_closed_r, accuracy


def find_last_closed_field_lines(kamodo_object, time_hours,
                                n_traces=36, r_start=3.0, r_max=100.0, r_step=8.0,
                                tolerance=0.1, coord_system='GSM',
                                use_parallel=False, use_numba=True, n_jobs=-1,
                                use_batches=None, max_memory_mb=2000,
                                show_progress=True, verbose=False):
    """
    Find last closed magnetic field lines and approximate magnetopause boundary.
    
    Uses optimized batch processing and parallel computation for improved performance.
    
    Parameters:
    -----------
    kamodo_object : object
        Kamodo object containing magnetic field data
    time_hours : float
        Time in hours since midnight of first day
    n_traces : int, optional (default=36)
        Number of longitude traces (evenly spaced around equator)
    r_start : float, optional (default=3.0)
        Starting radius for search in RE
    r_max : float, optional (default=100.0)
        Maximum radius to search in RE
    r_step : float, optional (default=8.0)
        Initial step size for coarse search in RE
    tolerance : float, optional (default=0.1)
        Tolerance for bisection method in RE
    coord_system : str, optional (default='GSM')
        Coordinate system for field tracing
    use_parallel : bool, optional (default=False)
        Whether to use parallel processing
    use_numba : bool, optional (default=True)
        Whether to use numba optimization
    n_jobs : int, optional (default=-1)
        Number of parallel jobs (-1 uses all cores)
    use_batches : bool, optional
        Whether to use batch processing (auto-determined if None)
    max_memory_mb : float, optional (default=2000)
        Maximum memory usage for batch processing in MB
    show_progress : bool, optional (default=True)
        Whether to show progress bars
    verbose : bool, optional (default=False)
        Whether to print detailed progress information

    Returns:
    --------
    results : dict
        Dictionary containing:
        - 'closed_traces': List of closed field line traces
        - 'open_traces': List of open field line traces
        - 'last_closed_r': Array of last closed radial distances for each longitude
        - 'longitudes': Array of longitude values in degrees
        - 'magnetopause_points': Array of approximate magnetopause boundary points
        - 'boundary_accuracy': Array of estimated accuracy for each boundary point
        - 'performance_stats': Dictionary with timing and performance information
        - 'summary': Dictionary with summary statistics
    """
    if verbose:
        print(f"Finding last closed field lines at t={time_hours:.2f} hours")
        print(f"Using bisection method with tolerance {tolerance:.3f} RE")
        print(f"Searching from {r_start} to {r_max} RE")
        print(f"Using {n_traces} traces at {360/n_traces:.1f}° intervals")
        print(f"Optimization settings: parallel={use_parallel}, numba={use_numba}")

    start_time = time.time()

    # Create optimized magnetic field tracer
    mag_tracer = vft.create_magnetic_tracer(
        kamodo_object,
        use_numba=use_numba,
        use_parallel=use_parallel,
        n_jobs=n_jobs,
        verbose=verbose
    )

    # Generate longitude angles (evenly spaced)
    longitudes = np.linspace(0, 360, n_traces, endpoint=False)

    # Initialize results storage
    closed_traces = []
    open_traces = []
    last_closed_r = np.full(n_traces, np.nan)
    boundary_accuracy = np.full(n_traces, np.nan)
    magnetopause_points = []

    if verbose:
        progress_interval = max(1, n_traces // 10)

    # Process each longitude
    print(f"Processing {n_traces} longitude slices...")
    
    # Use tqdm for progress if available and requested
    longitude_iter = enumerate(longitudes)
    if show_progress and vft.TQDM_AVAILABLE:
        from tqdm import tqdm
        longitude_iter = tqdm(longitude_iter, total=n_traces, desc="Finding last closed lines")

    for i, longitude_deg in longitude_iter:
        if verbose and (i % progress_interval == 0):
            print(f"Processing longitude {longitude_deg:.1f}° ({i+1}/{n_traces})")

        # Find last closed radius at this longitude using bisection
        last_r, accuracy = _find_last_closed_radius_bisection(
            mag_tracer, longitude_deg, time_hours,
            r_min=r_start, r_max=r_max, tolerance=tolerance,
            verbose=verbose
        )
        
        last_closed_r[i] = last_r
        boundary_accuracy[i] = accuracy

        # Create magnetopause boundary point
        if not np.isnan(last_r):
            longitude_rad = np.deg2rad(longitude_deg)
            boundary_point = [
                last_r * np.cos(longitude_rad),
                last_r * np.sin(longitude_rad),
                0.0
            ]
            magnetopause_points.append(boundary_point)

    # Convert magnetopause points to array
    magnetopause_points = np.array(magnetopause_points) if magnetopause_points else np.array([]).reshape(0, 3)

    # Now trace some example field lines for visualization
    if verbose:
        print("\nTracing example field lines for visualization...")

    # Select a subset of positions for detailed tracing
    n_example_traces = min(18, n_traces)  # Limit for visualization
    example_indices = np.linspace(0, n_traces-1, n_example_traces, dtype=int)
    
    # Prepare start points for batch tracing
    closed_start_points = []
    open_start_points = []
    
    for idx in example_indices:
        longitude_deg = longitudes[idx]
        longitude_rad = np.deg2rad(longitude_deg)
        last_r = last_closed_r[idx]
        
        if not np.isnan(last_r):
            # Add a closed field line (slightly inside last closed)
            r_closed = max(r_start, last_r - 0.5)
            closed_start_points.append([
                r_closed * np.cos(longitude_rad),
                r_closed * np.sin(longitude_rad),
                0.0
            ])
            
            # Add an open field line (slightly outside last closed)
            r_open = last_r + 0.5
            open_start_points.append([
                r_open * np.cos(longitude_rad),
                r_open * np.sin(longitude_rad),
                0.0
            ])

    # Batch trace closed field lines
    if closed_start_points:
        if verbose:
            print(f"Batch tracing {len(closed_start_points)} closed field lines...")
        
        closed_traces = mag_tracer.trace_multiple_lines(
            closed_start_points,
            time_hours,
            max_steps=3000,
            step_size_re=0.1,
            show_progress=show_progress and not vft.TQDM_AVAILABLE,
            use_batches=use_batches,
            max_memory_mb=max_memory_mb
        )

    # Batch trace open field lines
    if open_start_points:
        if verbose:
            print(f"Batch tracing {len(open_start_points)} open field lines...")
        
        open_traces = mag_tracer.trace_multiple_lines(
            open_start_points,
            time_hours,
            max_steps=3000,
            step_size_re=0.1,
            show_progress=show_progress and not vft.TQDM_AVAILABLE,
            use_batches=use_batches,
            max_memory_mb=max_memory_mb
        )

    # Calculate execution time and performance stats
    total_time = time.time() - start_time
    tracer_stats = mag_tracer.get_performance_stats()
    
    performance_stats = {
        'total_execution_time': total_time,
        'n_longitudes_processed': n_traces,
        'n_closed_traces': len(closed_traces),
        'n_open_traces': len(open_traces),
        'longitudes_per_second': n_traces / total_time,
        'tracer_stats': tracer_stats
    }

    # Create summary statistics
    valid_radii = last_closed_r[~np.isnan(last_closed_r)]
    summary = {
        'n_longitudes': n_traces,
        'n_valid_boundaries': len(valid_radii),
        'mean_boundary_radius': np.mean(valid_radii) if len(valid_radii) > 0 else np.nan,
        'std_boundary_radius': np.std(valid_radii) if len(valid_radii) > 0 else np.nan,
        'min_boundary_radius': np.min(valid_radii) if len(valid_radii) > 0 else np.nan,
        'max_boundary_radius': np.max(valid_radii) if len(valid_radii) > 0 else np.nan,
        'mean_accuracy': np.mean(boundary_accuracy[~np.isnan(boundary_accuracy)]),
        'execution_time': total_time
    }

    if verbose:
        print(f"\nCompleted in {total_time:.2f} seconds")
        print(f"Found {len(valid_radii)}/{n_traces} valid boundary points")
        if len(valid_radii) > 0:
            print(f"Mean boundary radius: {summary['mean_boundary_radius']:.2f} ± {summary['std_boundary_radius']:.2f} RE")
            print(f"Range: {summary['min_boundary_radius']:.2f} - {summary['max_boundary_radius']:.2f} RE")

    return {
        'closed_traces': closed_traces,
        'open_traces': open_traces,
        'last_closed_r': last_closed_r,
        'longitudes': longitudes,
        'magnetopause_points': magnetopause_points,
        'boundary_accuracy': boundary_accuracy,
        'performance_stats': performance_stats,
        'summary': summary
    }


def plot_last_closed_field_lines(results, backend='matplotlib', show_earth=True, 
                                interactive=True, figsize=(16, 12)):
    """
    Plot last closed field lines and magnetopause boundary.
    
    Parameters:
    -----------
    results : dict
        Results from find_last_closed_field_lines
    backend : str, optional (default='matplotlib')
        Plotting backend ('matplotlib' or 'plotly')
    show_earth : bool, optional (default=True)
        Whether to show Earth
    interactive : bool, optional (default=True)
        Whether to make plots interactive (for plotly)
    figsize : tuple, optional (default=(16, 12))
        Figure size for matplotlib
        
    Returns:
    --------
    fig : matplotlib.figure.Figure or plotly.graph_objects.Figure
        The generated figure
    """
    if backend.lower() == 'plotly' and PLOTLY_AVAILABLE:
        return _plot_last_closed_field_lines_plotly(results, show_earth, interactive)
    else:
        return _plot_last_closed_field_lines_matplotlib(results, show_earth, figsize)


def _plot_last_closed_field_lines_matplotlib(results, show_earth=True, figsize=(16, 12)):
    """
    Create matplotlib plots of last closed field lines.
    """
    fig = plt.figure(figsize=figsize)
    
    # Create subplots: 3D view, XY plane, boundary profile, performance stats
    ax1 = fig.add_subplot(221, projection='3d')  # 3D view
    ax2 = fig.add_subplot(222)  # XY plane
    ax3 = fig.add_subplot(223)  # Boundary radius vs longitude
    ax4 = fig.add_subplot(224)  # Performance stats

    # Extract data
    closed_traces = results['closed_traces']
    open_traces = results['open_traces']
    magnetopause_points = results['magnetopause_points']
    longitudes = results['longitudes']
    last_closed_r = results['last_closed_r']
    summary = results['summary']
    
    # Draw Earth
    if show_earth:
        r_earth = 1.0
        
        # Earth in 3D
        u = np.linspace(0, 2 * np.pi, 50)
        v = np.linspace(0, np.pi, 50)
        x = r_earth * np.outer(np.cos(u), np.sin(v))
        y = r_earth * np.outer(np.sin(u), np.sin(v))
        z = r_earth * np.outer(np.ones(np.size(u)), np.cos(v))
        ax1.plot_surface(x, y, z, color='lightblue', alpha=0.3)
        
        # Earth in XY plane
        circle = plt.Circle((0, 0), r_earth, color='lightblue', alpha=0.3)
        ax2.add_patch(circle)

    # Plot closed field lines in blue
    for trace in closed_traces:
        if 'combined' in trace:
            coords = trace['combined']['coordinates']
        elif 'forward' in trace:
            coords = trace['forward']['coordinates']
        else:
            continue
            
        ax1.plot(coords[:, 0], coords[:, 1], coords[:, 2], 'b-', alpha=0.7, linewidth=1)
        ax2.plot(coords[:, 0], coords[:, 1], 'b-', alpha=0.7, linewidth=1)

    # Plot open field lines in red
    for trace in open_traces:
        if 'combined' in trace:
            coords = trace['combined']['coordinates']
        elif 'forward' in trace:
            coords = trace['forward']['coordinates']
        else:
            continue
            
        ax1.plot(coords[:, 0], coords[:, 1], coords[:, 2], 'r-', alpha=0.7, linewidth=1)
        ax2.plot(coords[:, 0], coords[:, 1], 'r-', alpha=0.7, linewidth=1)

    # Plot magnetopause boundary
    if len(magnetopause_points) > 0:
        # Close the boundary by connecting last point to first
        boundary_closed = np.vstack([magnetopause_points, magnetopause_points[0]])
        
        ax1.plot(boundary_closed[:, 0], boundary_closed[:, 1], boundary_closed[:, 2], 
                'g-', linewidth=3, label='Magnetopause')
        ax2.plot(boundary_closed[:, 0], boundary_closed[:, 1], 
                'g-', linewidth=3, label='Magnetopause')
        ax2.scatter(magnetopause_points[:, 0], magnetopause_points[:, 1], 
                   c='green', s=30, alpha=0.7)

    # Plot boundary radius vs longitude
    valid_mask = ~np.isnan(last_closed_r)
    if np.any(valid_mask):
        ax3.plot(longitudes[valid_mask], last_closed_r[valid_mask], 'go-', 
                linewidth=2, markersize=4)
        ax3.fill_between(longitudes[valid_mask], 
                        last_closed_r[valid_mask] - results['boundary_accuracy'][valid_mask],
                        last_closed_r[valid_mask] + results['boundary_accuracy'][valid_mask],
                        alpha=0.3, color='green')

    # Performance statistics
    perf_stats = results['performance_stats']
    stats_text = [
        f"Total time: {perf_stats['total_execution_time']:.2f} s",
        f"Longitudes/sec: {perf_stats['longitudes_per_second']:.2f}",
        f"Closed traces: {perf_stats['n_closed_traces']}",
        f"Open traces: {perf_stats['n_open_traces']}",
        f"Valid boundaries: {summary['n_valid_boundaries']}/{summary['n_longitudes']}",
        f"Mean radius: {summary['mean_boundary_radius']:.2f} ± {summary['std_boundary_radius']:.2f} RE",
        f"Range: {summary['min_boundary_radius']:.2f} - {summary['max_boundary_radius']:.2f} RE"
    ]
    
    ax4.text(0.05, 0.95, '\n'.join(stats_text), transform=ax4.transAxes, 
            verticalalignment='top', fontfamily='monospace', fontsize=10)
    ax4.set_xlim(0, 1)
    ax4.set_ylim(0, 1)
    ax4.axis('off')

    # Set labels and titles
    ax1.set_xlabel('X (RE)')
    ax1.set_ylabel('Y (RE)')
    ax1.set_zlabel('Z (RE)')
    ax1.set_title('3D View: Last Closed Magnetic Field Lines')

    ax2.set_xlabel('X (RE)')
    ax2.set_ylabel('Y (RE)')
    ax2.set_title('XY Plane View')
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    ax3.set_xlabel('Longitude (degrees)')
    ax3.set_ylabel('Last Closed Radius (RE)')
    ax3.set_title('Magnetopause Boundary Profile')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, 360)

    ax4.set_title('Performance Statistics')

    plt.tight_layout()
    return fig


def _plot_last_closed_field_lines_plotly(results, show_earth=True, interactive=True):
    """
    Create plotly plots of last closed field lines.
    """
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{'type': 'scene'}, {'type': 'xy'}],
               [{'type': 'xy'}, {'type': 'xy'}]],
        subplot_titles=('3D View: Last Closed Magnetic Field Lines',
                       'XY Plane View',
                       'Magnetopause Boundary Profile',
                       'Performance Summary')
    )

    # Extract data
    closed_traces = results['closed_traces']
    open_traces = results['open_traces']
    magnetopause_points = results['magnetopause_points']
    longitudes = results['longitudes']
    last_closed_r = results['last_closed_r']
    summary = results['summary']

    # Draw Earth
    if show_earth:
        r_earth = 1.0
        
        # Earth sphere in 3D
        theta = np.linspace(0, 2 * np.pi, 50)
        phi = np.linspace(0, np.pi, 50)
        x = r_earth * np.outer(np.cos(theta), np.sin(phi))
        y = r_earth * np.outer(np.sin(theta), np.sin(phi))
        z = r_earth * np.outer(np.ones(np.size(theta)), np.cos(phi))
        
        fig.add_trace(
            go.Surface(x=x, y=y, z=z, 
                      colorscale=[[0, 'lightblue'], [1, 'lightblue']],
                      opacity=0.3, showscale=False, name='Earth'),
            row=1, col=1
        )
        
        # Earth circle in XY plane
        theta_circle = np.linspace(0, 2*np.pi, 100)
        x_circle = r_earth * np.cos(theta_circle)
        y_circle = r_earth * np.sin(theta_circle)
        
        fig.add_trace(
            go.Scatter(x=x_circle, y=y_circle, mode='lines',
                      line=dict(color='lightblue', width=2),
                      fill='toself', fillcolor='rgba(173, 216, 230, 0.3)',
                      showlegend=False, name='Earth'),
            row=1, col=2
        )

    # Plot closed field lines
    for i, trace in enumerate(closed_traces):
        if 'combined' in trace:
            coords = trace['combined']['coordinates']
        elif 'forward' in trace:
            coords = trace['forward']['coordinates']
        else:
            continue
        
        # 3D plot
        fig.add_trace(
            go.Scatter3d(x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
                        mode='lines', line=dict(color='blue', width=2),
                        showlegend=(i==0), name='Closed field lines'),
            row=1, col=1
        )
        
        # XY plane
        fig.add_trace(
            go.Scatter(x=coords[:, 0], y=coords[:, 1], mode='lines',
                      line=dict(color='blue', width=2),
                      showlegend=False),
            row=1, col=2
        )

    # Plot open field lines
    for i, trace in enumerate(open_traces):
        if 'combined' in trace:
            coords = trace['combined']['coordinates']
        elif 'forward' in trace:
            coords = trace['forward']['coordinates']
        else:
            continue
        
        # 3D plot
        fig.add_trace(
            go.Scatter3d(x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
                        mode='lines', line=dict(color='red', width=2),
                        showlegend=(i==0), name='Open field lines'),
            row=1, col=1
        )
        
        # XY plane
        fig.add_trace(
            go.Scatter(x=coords[:, 0], y=coords[:, 1], mode='lines',
                      line=dict(color='red', width=2),
                      showlegend=False),
            row=1, col=2
        )

    # Plot magnetopause boundary
    if len(magnetopause_points) > 0:
        # Close the boundary
        boundary_closed = np.vstack([magnetopause_points, magnetopause_points[0]])
        
        # 3D boundary
        fig.add_trace(
            go.Scatter3d(x=boundary_closed[:, 0], y=boundary_closed[:, 1], z=boundary_closed[:, 2],
                        mode='lines', line=dict(color='green', width=4),
                        name='Magnetopause'),
            row=1, col=1
        )
        
        # XY boundary
        fig.add_trace(
            go.Scatter(x=boundary_closed[:, 0], y=boundary_closed[:, 1],
                      mode='lines+markers', 
                      line=dict(color='green', width=4),
                      marker=dict(color='green', size=6),
                      showlegend=False),
            row=1, col=2
        )

    # Boundary profile
    valid_mask = ~np.isnan(last_closed_r)
    if np.any(valid_mask):
        # Main boundary line
        fig.add_trace(
            go.Scatter(x=longitudes[valid_mask], y=last_closed_r[valid_mask],
                      mode='lines+markers',
                      line=dict(color='green', width=3),
                      marker=dict(color='green', size=6),
                      name='Boundary radius',
                      showlegend=False),
            row=2, col=1
        )
        
        # Error band
        upper_bound = last_closed_r[valid_mask] + results['boundary_accuracy'][valid_mask]
        lower_bound = last_closed_r[valid_mask] - results['boundary_accuracy'][valid_mask]
        
        fig.add_trace(
            go.Scatter(x=np.concatenate([longitudes[valid_mask], longitudes[valid_mask][::-1]]),
                      y=np.concatenate([upper_bound, lower_bound[::-1]]),
                      fill='toself', fillcolor='rgba(0,255,0,0.3)',
                      line=dict(color='rgba(255,255,255,0)'),
                      showlegend=False),
            row=2, col=1
        )

    # Performance summary
    perf_stats = results['performance_stats']
    summary_text = (
        f"<b>Performance Summary</b><br><br>"
        f"Total execution time: {perf_stats['total_execution_time']:.2f} s<br>"
        f"Longitudes processed/sec: {perf_stats['longitudes_per_second']:.2f}<br>"
        f"Closed field lines traced: {perf_stats['n_closed_traces']}<br>"
        f"Open field lines traced: {perf_stats['n_open_traces']}<br><br>"
        f"<b>Boundary Statistics</b><br><br>"
        f"Valid boundary points: {summary['n_valid_boundaries']}/{summary['n_longitudes']}<br>"
        f"Mean radius: {summary['mean_boundary_radius']:.2f} ± {summary['std_boundary_radius']:.2f} RE<br>"
        f"Range: {summary['min_boundary_radius']:.2f} - {summary['max_boundary_radius']:.2f} RE<br>"
        f"Mean accuracy: {summary['mean_accuracy']:.3f} RE"
    )
    
    fig.add_annotation(
        text=summary_text,
        xref="paper", yref="paper", xshift=0, yshift=100,
        x=0.5, y=1.0, xanchor='center', yanchor='top',
        showarrow=False,
        font=dict(size=12, family="monospace"),
        row=2, col=2
    )

    # Update layout
    fig.update_scenes(
        xaxis_title='X (RE)', yaxis_title='Y (RE)', zaxis_title='Z (RE)',
        aspectmode='data',
        row=1, col=1
    )
    
    fig.update_xaxes(title_text='X (RE)', scaleanchor='y', scaleratio=1, row=1, col=2)
    fig.update_yaxes(title_text='Y (RE)', row=1, col=2)
    
    fig.update_xaxes(title_text='Longitude (degrees)', range=[0, 360], row=2, col=1)
    fig.update_yaxes(title_text='Last Closed Radius (RE)', row=2, col=1)
    
    # Empty axes for the summary panel
    fig.update_xaxes(showticklabels=False, showgrid=False, zeroline=False, row=2, col=2)
    fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False, row=2, col=2)
    
    # Overall layout settings
    fig.update_layout(
        title='Last Closed Magnetic Field Lines and Magnetopause Boundary',
        height=800,
        width=1200,
        hovermode='closest' if interactive else False,
        legend=dict(
            x=0.01, y=0.99,
            bgcolor='rgba(255,255,255,0.7)',
            bordercolor='rgba(0,0,0,0.5)',
            borderwidth=1
        )
    )
    
    return fig


def find_last_closed_field_lines_auto_bounds(ko, time_hours, n_traces=36,
                                            r_start=2.0, r_max_factor=0.8,
                                            tolerance=0.01, max_bisection_iter=10,
                                            initial_step=1.0, max_steps=5000,
                                            step_size_re=0.1, coord_system='GSM',
                                            use_parallel=True, use_numba=True, n_jobs=-1,
                                            use_batches=None, max_memory_mb=2000,
                                            show_progress=True, verbose=False):
    """
    Find last closed magnetic field lines with automatic bound detection.
    
    Parameters:
    -----------
    ko : Kamodo object
        Kamodo object containing magnetic field data
    time_hours : float
        Time in hours since midnight of first day of data
    n_traces : int
        Number of traces (default 36 for 10-degree spacing)
    r_start : float
        Starting radius, minimum 2.0 RE (default 2.0 RE)
    r_max_factor : float
        Maximum search radius as fraction of model extent (default 0.8)
    tolerance : float
        Tolerance for bisection method in RE (default 0.01 RE)
    max_bisection_iter : int
        Maximum iterations for bisection method (default 10)
    initial_step : float
        Initial step size for coarse search in RE (default 1.0 RE)
    max_steps : int
        Maximum integration steps per trace
    step_size_re : float
        Integration step size in RE
    coord_system : str
        Coordinate system ('GSM', 'GSE', 'SM', etc.)
    use_parallel : bool
        Whether to use parallel processing
    use_numba : bool
        Whether to use numba optimization
    n_jobs : int
        Number of parallel jobs (-1 uses all cores)
    use_batches : bool, optional
        Whether to use batch processing (auto-determined if None)
    max_memory_mb : float
        Maximum memory usage for batch processing in MB
    show_progress : bool
        Whether to show progress bars
    verbose : bool
        Whether to print detailed information during execution
    
    Returns:
    --------
    results : dict
        Dictionary containing found field lines and boundary points
    """
    # Create tracer to detect model bounds
    test_tracer = vft.create_magnetic_tracer(
        ko,
        use_numba=use_numba,
        use_parallel=use_parallel,
        n_jobs=n_jobs,
        verbose=False
    )

    # Try to determine model bounds by sampling test points
    try:
        bounds = {'x_range': [-30, 30], 'y_range': [-30, 30], 'z_range': [-30, 30]}
        test_radius = 25.0
        test_points = []
        
        # Generate test points at different radii
        for r in [5.0, 10.0, 15.0, 20.0, test_radius]:
            for theta in [0, np.pi/2, np.pi, 3*np.pi/2]:
                test_points.append([r * np.cos(theta), r * np.sin(theta), 0.0])
        
        # Test which points are within model bounds
        valid_points = []
        for point in test_points:
            try:
                test_trace = test_tracer.trace_vector_line(
                    point, time_hours, max_steps=10
                )
                if test_trace is not None:
                    valid_points.append(point)
            except Exception:
                pass  # If error, assume point is outside model bounds
                
        # Estimate bounds from valid points
        if valid_points:
            valid_points = np.array(valid_points)
            max_r = np.max(np.sqrt(np.sum(valid_points[:, :2]**2, axis=1)))
            bounds = {
                'x_range': [-max_r, max_r],
                'y_range': [-max_r, max_r],
                'z_range': [-10, 10]  # Conservative z-range
            }
            if verbose:
                print(f"Detected model bounds: r_max={max_r:.1f} RE")
        else:
            if verbose:
                print(f"Could not detect model bounds, using defaults")
    except Exception as e:
        if verbose:
            print(f"Error detecting bounds: {e}, using defaults")
    
    # Calculate maximum radius for searching
    x_extent = bounds['x_range'][1] - bounds['x_range'][0]
    y_extent = bounds['y_range'][1] - bounds['y_range'][0]
    model_radius = min(x_extent, y_extent) / 2
    r_max = min(r_max_factor * model_radius, 20.0)  # Limit to 20 RE
    
    if verbose:
        print(f"Using r_max = {r_max:.1f} RE (model_radius = {model_radius:.1f} RE)")
    
    # Call the standard function with the determined bounds
    results = find_last_closed_field_lines(
        ko, time_hours,
        n_traces=n_traces,
        r_start=r_start,
        r_max=r_max,
        r_step=initial_step,
        tolerance=tolerance,
        coord_system=coord_system,
        use_parallel=use_parallel,
        use_numba=use_numba,
        n_jobs=n_jobs,
        use_batches=use_batches,
        max_memory_mb=max_memory_mb,
        show_progress=show_progress,
        verbose=verbose
    )
    
    # Add bounds information
    results['model_bounds'] = bounds
    results['search_parameters'] = {
        'r_start': r_start,
        'r_max': r_max,
        'r_max_factor': r_max_factor
    }
    
    return results


def analyze_magnetopause_shape(results):
    """
    Analyze the shape of the magnetopause boundary.
    
    Parameters:
    -----------
    results : dict
        Results from find_last_closed_field_lines
        
    Returns:
    --------
    analysis : dict
        Dictionary with magnetopause shape analysis
    """
    magnetopause_points = results['magnetopause_points']
    longitudes = results['longitudes']
    last_closed_r = results['last_closed_r']
    
    if len(magnetopause_points) == 0 or not np.any(~np.isnan(last_closed_r)):
        return {
            'status': 'error',
            'message': 'No valid magnetopause boundary points found'
        }
    
    # Filter valid points
    valid_mask = ~np.isnan(last_closed_r)
    valid_longitudes = longitudes[valid_mask]
    valid_radii = last_closed_r[valid_mask]
    valid_points = magnetopause_points
    
    # 1. Calculate basic statistics
    mean_radius = np.mean(valid_radii)
    std_radius = np.std(valid_radii)
    min_radius = np.min(valid_radii)
    min_long = valid_longitudes[np.argmin(valid_radii)]
    max_radius = np.max(valid_radii)
    max_long = valid_longitudes[np.argmax(valid_radii)]
    
    # 2. Day/night asymmetry
    day_mask = (valid_longitudes >= 270) | (valid_longitudes <= 90)
    night_mask = ~day_mask
    
    day_mean = np.mean(valid_radii[day_mask]) if np.any(day_mask) else np.nan
    night_mean = np.mean(valid_radii[night_mask]) if np.any(night_mask) else np.nan
    day_night_ratio = day_mean / night_mean if not np.isnan(night_mean) and night_mean > 0 else np.nan
    
    # 3. Dawn/dusk asymmetry
    dawn_mask = (valid_longitudes >= 0) & (valid_longitudes <= 180)
    dusk_mask = ~dawn_mask
    
    dawn_mean = np.mean(valid_radii[dawn_mask]) if np.any(dawn_mask) else np.nan
    dusk_mean = np.mean(valid_radii[dusk_mask]) if np.any(dusk_mask) else np.nan
    dawn_dusk_ratio = dawn_mean / dusk_mean if not np.isnan(dusk_mean) and dusk_mean > 0 else np.nan
    
    # 4. Subsolar point (approximation)
    subsolar_idx = np.abs(valid_longitudes - 0).argmin()
    subsolar_r = valid_radii[subsolar_idx]
    
    return {
        'status': 'success',
        'n_boundary_points': len(valid_points),
        'mean_radius': mean_radius,
        'std_radius': std_radius,
        'min_radius': {'value': min_radius, 'longitude': min_long},
        'max_radius': {'value': max_radius, 'longitude': max_long},
        'subsolar_point': {'radius': subsolar_r, 'longitude': valid_longitudes[subsolar_idx]},
        'asymmetry': {
            'day_mean': day_mean,
            'night_mean': night_mean,
            'day_night_ratio': day_night_ratio,
            'dawn_mean': dawn_mean,
            'dusk_mean': dusk_mean,
            'dawn_dusk_ratio': dawn_dusk_ratio,
        }
    }


def benchmark_last_closed_finder(ko, time_hours, n_traces=18,
                               use_numba_options=[True, False], 
                               use_parallel_options=[True, False],
                               verbose=True):
    """
    Benchmark different configurations of the last closed field line finder.
    
    Parameters:
    -----------
    ko : Kamodo object
        Kamodo object with magnetic field components
    time_hours : float
        Time in hours
    n_traces : int
        Number of longitude slices to trace
    use_numba_options : list
        List of boolean options to test for numba usage
    use_parallel_options : list
        List of boolean options to test for parallel processing
    verbose : bool
        Whether to print detailed information
        
    Returns:
    --------
    results : dict
        Dictionary of benchmark results
    """
    print(f"Starting benchmark with {n_traces} longitude slices")
    
    results = {}
    
    for use_numba in use_numba_options:
        for use_parallel in use_parallel_options:
            # Skip if required features aren't available
            if use_numba and not vft.NUMBA_AVAILABLE:
                continue
            if use_parallel and not vft.JOBLIB_AVAILABLE:
                continue
                
            config_name = f"numba={use_numba}_parallel={use_parallel}"
            print(f"\nTesting configuration: {config_name}")
            
            start_time = time.time()
            
            # Run the finder with limited scope for benchmark
            result = find_last_closed_field_lines(
                ko, time_hours,
                n_traces=n_traces,
                r_start=2.0,
                r_max=12.0,
                r_step=1.0,
                tolerance=0.2,  # Use coarse tolerance for speed
                use_numba=use_numba,
                use_parallel=use_parallel,
                verbose=verbose
            )
            
            end_time = time.time()
            execution_time = end_time - start_time
            
            # Store results
            results[config_name] = {
                'execution_time': execution_time,
                'valid_boundaries': result['summary']['n_valid_boundaries'],
                'mean_radius': result['summary']['mean_boundary_radius'],
                'longitudes_per_second': n_traces / execution_time
            }
            
            print(f"Execution time: {execution_time:.2f} seconds")
            print(f"Performance: {n_traces / execution_time:.2f} longitudes/second")
            print(f"Valid boundaries: {result['summary']['n_valid_boundaries']}/{n_traces}")
    
    # Print summary
    if results:
        print("\n===== BENCHMARK SUMMARY =====")
        baseline = results.get("numba=False_parallel=False", None)
        baseline_time = baseline['execution_time'] if baseline else None
        
        for config_name, data in results.items():
            speedup = baseline_time / data['execution_time'] if baseline_time else 1.0
            print(f"{config_name}: {data['execution_time']:.2f} sec, "
                  f"{speedup:.2f}x speedup, "
                  f"{data['longitudes_per_second']:.2f} longitudes/sec")
    
    return results


def example_find_last_closed_field_lines(ko, time_hours=12.0):
    """
    Example of finding and plotting last closed field lines.
    
    Parameters:
    -----------
    ko : Kamodo object
        Your Kamodo object
    time_hours : float
        Time in hours since midnight of first day
    """
    print("Example: Finding Last Closed Field Lines")
    print("=" * 40)
    
    # Find last closed field lines using auto-bounds
    results = find_last_closed_field_lines_auto_bounds(
        ko, time_hours,
        n_traces=36,      # 36 traces at 10° intervals
        r_start=2.0,      # Start searching at 2 RE
        verbose=True,
        show_progress=True,
        use_numba=True,
        use_parallel=True
    )
    
    # Analyze the shape
    shape_analysis = analyze_magnetopause_shape(results)
    
    # Create plots
    print("\nCreating plots...")
    
    # Matplotlib version
    fig_mpl = plot_last_closed_field_lines(results, backend='matplotlib')
    plt.show()
    
    # Plotly version (interactive) - if available
    if PLOTLY_AVAILABLE:
        fig_plotly = plot_last_closed_field_lines(results, backend='plotly')
        
        # Show plotly in notebook if in notebook environment
        try:
            from IPython import get_ipython
            if get_ipython() is not None:
                fig_plotly.show()
        except ImportError:
            print("To view interactive plot, call: fig_plotly.show()")
    
    # Print summary
    summary = results['summary']
    print(f"\nResults Summary:")
    print(f"  Execution time: {summary['execution_time']:.2f} seconds")
    print(f"  Valid boundary points: {summary['n_valid_boundaries']}/{summary['n_longitudes']}")
    print(f"  Mean boundary radius: {summary['mean_boundary_radius']:.2f} ± {summary['std_boundary_radius']:.2f} RE")
    print(f"  Range: {summary['min_boundary_radius']:.2f} - {summary['max_boundary_radius']:.2f} RE")
    
    if shape_analysis['status'] == 'success':
        print("\nMagnetopause Shape Analysis:")
        asym = shape_analysis['asymmetry']
        print(f"  Subsolar point: {shape_analysis['subsolar_point']['radius']:.2f} RE")
        print(f"  Day/night asymmetry: {asym['day_night_ratio']:.2f}")
        print(f"  Dawn/dusk asymmetry: {asym['dawn_dusk_ratio']:.2f}")
    
    return results, fig_mpl


if __name__ == "__main__":
    # This code runs when the script is executed directly
    print("Kamodo Last Closed Field Line Finder")
    print("To use this module, import it and call its functions.")
    print("For an example, use:")
    print("  results = find_last_closed_field_lines(kamodo_object, time_hours)")
    print("  fig = plot_last_closed_field_lines(results)")
    
    # Check available optimizations
    print("\nAvailable optimizations:")
    print(f"  Numba JIT compilation: {'✓' if vft.NUMBA_AVAILABLE else '✗'}")
    print(f"  Parallel processing:   {'✓' if vft.JOBLIB_AVAILABLE else '✗'}")
    print(f"  Progress bars (tqdm):  {'✓' if vft.TQDM_AVAILABLE else '✗'}")
    
    if not vft.NUMBA_AVAILABLE:
        print("\n⚠️  For best performance, install numba:")
        print("    pip install numba")
    
    if not vft.JOBLIB_AVAILABLE:
        print("\n⚠️  For parallel processing, install joblib:")
        print("    pip install joblib")

