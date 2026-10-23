"""
Kamodo Vector Field Tracer - Core field line tracing functionality

High-performance vector field tracer for Kamodo models.
Supports magnetic field lines, velocity streamlines, current streamlines, etc.
All coordinates are in Earth radii (RE).
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import gc
import kamodo_ccmc.flythrough.model_wrapper as MW

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    def tqdm(iterable=None, **kwargs):
        return iterable
    TQDM_AVAILABLE = False

def get_colormap(name):
    try:
        import matplotlib
        if hasattr(matplotlib, 'colormaps'):
            return matplotlib.colormaps[name]
        else:
            return plt.cm.get_cmap(name)
    except (AttributeError, KeyError):
        try:
            return plt.get_cmap(name)
        except:
            return plt.get_cmap('viridis')

class KamodoVectorFieldTracer:
    def __init__(self, kamodo_object, field_type='magnetic', component_names=None, 
                 verbose=False, zeroX=False, zeroY=False, zeroZ=False):
        self.kamodo_object = kamodo_object
        self.field_type = field_type.lower()
        self.verbose = verbose
        self.xyztype = 'xyz'
        self.zeroX = zeroX
        self.zeroY = zeroY
        self.zeroZ = zeroZ
        
        if component_names is None:
            prefix = {'magnetic': 'B', 'velocity': 'v', 'current': 'J', 'electric': 'E'}.get(self.field_type, 'f')
            self.component_names = [f'{prefix}_x', f'{prefix}_y', f'{prefix}_z']
        else:
            self.component_names = component_names
        
        self._cache_field_components()
        
        self.default_params = {
            'max_steps': 5000,
            'step_size_re': 0.05,
            'direction': 'both',
            'min_altitude_re': 1.0,
            'max_altitude_re': 50.0,
            'adaptive_stepping': True,
            'batch_size': 1000,
            'solver': 'rk4'  # 'euler' matches original baseline, 'rk4' for accuracy
        }
        
        self.performance_stats = {
            'total_traces': 0,
            'total_time': 0.0,
            'total_points': 0
        }
        
        self._determine_call_method()
    
    def _cache_field_components(self):
        self.field_components = []
        for component in self.component_names:
            if component in self.kamodo_object:
                self.field_components.append(self.kamodo_object[component])
            else:
                raise ValueError(f"Component '{component}' not found in Kamodo object.")
    
    def _determine_call_method(self):
        self._call_method = 'array_with_time'
        t_c = self.component_names[0]
        t_cr = MW.Coord_Range(self.kamodo_object, [t_c], print_output=False, return_dict=True)
        if 'X' in t_cr[t_c]:
            self.xyztype = 'XYZ'
            
        if self.xyztype == 'xyz':
            coords = [t_cr[t_c]['time'][0], t_cr[t_c]['x'][0], t_cr[t_c]['y'][0], t_cr[t_c]['z'][0]]
        else:
            coords = [t_cr[t_c]['time'][0], t_cr[t_c]['X'][0], t_cr[t_c]['Y'][0], t_cr[t_c]['Z'][0]]
            
        try:
            value = np.array(self.kamodo_object[t_c](coords))
            if len(value.shape) == 1:
                return
        except Exception:
            pass
            
        self._call_method = 'gridded'
        try:
            if self.xyztype == 'xyz':
                value = np.array(self.kamodo_object[t_c](time=t_cr[t_c]['time'][0], x=t_cr[t_c]['x'][0], y=t_cr[t_c]['y'][0], z=t_cr[t_c]['z'][0]))
            else:
                value = np.array(self.kamodo_object[t_c](time=t_cr[t_c]['time'][0], x=t_cr[t_c]['X'][0], y=t_cr[t_c]['Y'][0], z=t_cr[t_c]['Z'][0]))
            if len(value.shape) == 0:
                return
        except Exception:
            pass
            
        self._call_method = 'array_with_time'

    def get_field_vectors_batch(self, positions_re, time_hours):
        """
        Get field vectors for multiple positions simultaneously using C-level Kamodo batching.
        """
        positions_re = np.asarray(positions_re, dtype=np.float64)
        if positions_re.ndim == 1:
            positions_re = positions_re.reshape(1, -1)
            
        n_points = len(positions_re)
        
        if self._call_method == 'gridded':
            ax, ay, az = positions_re[:, 0], positions_re[:, 1], positions_re[:, 2]
        else:
            t_arr = np.full((n_points, 1), time_hours)
            coords_array = np.hstack((t_arr, positions_re)).tolist()
        
        field_components_values = []
        for i, component in enumerate(self.field_components):
            if self._call_method == 'array_with_time':
                values = component(coords_array)
            else:
                values = component(time=time_hours, x=ax, y=ay, z=az)
            
            field_components_values.append(np.asarray(values, dtype=np.float64).flatten())
            
        field_vectors = np.column_stack(field_components_values)
        field_magnitudes = np.linalg.norm(field_vectors, axis=1)
        
        if self.zeroX: field_vectors[:, 0] = 0.0
        if self.zeroY: field_vectors[:, 1] = 0.0
        if self.zeroZ: field_vectors[:, 2] = 0.0
        
        return field_vectors, field_magnitudes

    def _trace_vectorized_batch(self, start_points_re, time_hours, forward=True):
        """
        Zero-loop vectorized field line tracing with adaptive stepping and memory masking.
        """
        n_lines = len(start_points_re)
        sign = 1.0 if forward else -1.0
        base_step = self.params['step_size_re']
        max_steps = self.params['max_steps']
        use_adaptive = self.params['adaptive_stepping']
        solver = self.params['solver']
        
        all_coords = np.zeros((n_lines, max_steps + 1, 3))
        all_coords[:, 0, :] = start_points_re
        
        all_mags = np.zeros((n_lines, max_steps + 1))
        all_comps = np.zeros((n_lines, max_steps + 1, 3))
        
        initial_fields, initial_mags = self.get_field_vectors_batch(start_points_re, time_hours)
        all_mags[:, 0] = initial_mags
        all_comps[:, 0, :] = initial_fields
        
        active_mask = np.ones(n_lines, dtype=bool)
        step_counts = np.zeros(n_lines, dtype=int)
        stop_reasons = np.full(n_lines, "max_steps_reached", dtype=object)
        current_positions = np.array(start_points_re, dtype=np.float64)
        
        steps = 0
        while np.any(active_mask) and steps < max_steps:
            active_positions = current_positions[active_mask]

            if use_adaptive:
                radii = np.linalg.norm(active_positions, axis=1)
                # Scale step size linearly with distance, capping at 1.5 RE in the far tail.
                # At R=10 RE, step is 0.5 RE. At R>=30 RE, step is 1.5 RE.
                # Minimum step is 0.05 RE to preserve dipole curvature near Earth.
                step_sizes = np.clip(0.05 * radii, 0.05, 1.5)[:, np.newaxis]
            else:
                step_sizes = base_step

            # K1 / Euler Eval
            fields_k1, mags_k1 = self.get_field_vectors_batch(active_positions, time_hours)
            safe_mags1 = np.maximum(mags_k1, 1e-10)[:, np.newaxis]
            k1 = sign * step_sizes * (fields_k1 / safe_mags1)
            
            if solver == 'euler':
                new_positions = active_positions + k1
            else:
                # Full RK4 
                pos_k2 = active_positions + 0.5 * k1
                fields_k2, mags_k2 = self.get_field_vectors_batch(pos_k2, time_hours)
                safe_mags2 = np.maximum(mags_k2, 1e-10)[:, np.newaxis]
                k2 = sign * step_sizes * (fields_k2 / safe_mags2)
                
                pos_k3 = active_positions + 0.5 * k2
                fields_k3, mags_k3 = self.get_field_vectors_batch(pos_k3, time_hours)
                safe_mags3 = np.maximum(mags_k3, 1e-10)[:, np.newaxis]
                k3 = sign * step_sizes * (fields_k3 / safe_mags3)
                
                pos_k4 = active_positions + k3
                fields_k4, mags_k4 = self.get_field_vectors_batch(pos_k4, time_hours)
                safe_mags4 = np.maximum(mags_k4, 1e-10)[:, np.newaxis]
                k4 = sign * step_sizes * (fields_k4 / safe_mags4)
                
                new_positions = active_positions + (k1 + 2*k2 + 2*k3 + k4) / 6.0
                
            all_coords[active_mask, steps + 1, :] = new_positions
            all_mags[active_mask, steps + 1] = mags_k1
            all_comps[active_mask, steps + 1, :] = fields_k1
            
            current_positions[active_mask] = new_positions
            step_counts[active_mask] += 1
            
            # C-Level Boundary Checks
            radii_new = np.linalg.norm(new_positions, axis=1)
            hit_min = radii_new < self.params['min_altitude_re']
            hit_max = radii_new > self.params['max_altitude_re']
            too_weak = mags_k1 < 1e-10
            
            new_stopped = hit_min | hit_max | too_weak
            if np.any(new_stopped):
                global_active = np.where(active_mask)[0]
                
                stop_reasons[global_active[hit_min]] = "hit_min_altitude"
                stop_reasons[global_active[hit_max]] = "hit_max_altitude"
                stop_reasons[global_active[too_weak & ~hit_min & ~hit_max]] = "field_too_weak"
                
                active_mask[global_active[new_stopped]] = False
            
            steps += 1
            
        results = []
        for i in range(n_lines):
            n_pts = step_counts[i] + 1
            coords_out = all_coords[i, :n_pts, :].copy()
            results.append({
                'coordinates': coords_out,
                'field_magnitude': all_mags[i, :n_pts].copy(),
                'field_components': all_comps[i, :n_pts, :].copy(),
                'length': np.sum(np.linalg.norm(np.diff(coords_out, axis=0), axis=1)) if n_pts > 1 else 0.0,
                'n_points': n_pts,
                'stop_reason': stop_reasons[i]
            })
            
        return results

    def trace_multiple_lines(self, start_points_re, time_hours, show_progress=False, use_batches=None, **kwargs):
        self.params = self.default_params.copy()
        self.params.update(kwargs)
        
        start_points_re = np.array(start_points_re)
        n_lines = len(start_points_re)
        
        if use_batches is None:
            use_batches = n_lines > self.params['batch_size']
            
        batch_size = self.params['batch_size'] if use_batches else n_lines
        n_batches = (n_lines + batch_size - 1) // batch_size

        if self.verbose:
            print(f"Vector tracing {n_lines} lines at t={time_hours:.2f} using {self.params['solver'].upper()} in {n_batches} batches...")

        trace_forward = self.params['direction'] in ['forward', 'both']
        trace_backward = self.params['direction'] in ['backward', 'both']
        
        all_results = []
        start_time = time.time()
        
        batch_iter = range(n_batches)
        if show_progress and TQDM_AVAILABLE:
            batch_iter = tqdm(batch_iter, desc="Processing Vector Batches", unit="batch")

        for i in batch_iter:
            start_idx = i * batch_size
            end_idx = min((i + 1) * batch_size, n_lines)
            batch_points = start_points_re[start_idx:end_idx]
            batch_n = len(batch_points)
            
            batch_results = [{} for _ in range(batch_n)]
            
            if trace_forward:
                forward_res = self._trace_vectorized_batch(batch_points, time_hours, forward=True)
                for j in range(batch_n): batch_results[j]['forward'] = forward_res[j]
                    
            if trace_backward:
                backward_res = self._trace_vectorized_batch(batch_points, time_hours, forward=False)
                for j in range(batch_n): batch_results[j]['backward'] = backward_res[j]
                    
            if trace_forward and trace_backward:
                for j in range(batch_n):
                    fwd = batch_results[j]['forward']
                    bwd = batch_results[j]['backward']
                    
                    rev_bwd_coords = bwd['coordinates'][::-1]
                    combined_coords = np.vstack((rev_bwd_coords[:-1], fwd['coordinates']))
                    
                    rev_bwd_mag = bwd['field_magnitude'][::-1]
                    combined_mag = np.concatenate((rev_bwd_mag[:-1], fwd['field_magnitude']))
                    
                    rev_bwd_comp = bwd['field_components'][::-1]
                    combined_comp = np.vstack((rev_bwd_comp[:-1], fwd['field_components']))
                    
                    batch_results[j]['combined'] = {
                        'coordinates': combined_coords,
                        'field_magnitude': combined_mag,
                        'field_components': combined_comp,
                        'length': bwd['length'] + fwd['length'],
                        'n_points': len(combined_coords),
                        'stop_reasons': [bwd['stop_reason'], fwd['stop_reason']]
                    }
            
            all_results.extend(batch_results)
            
            if i % 3 == 0 and i > 0: gc.collect()

        self.performance_stats['total_time'] += (time.time() - start_time)
        self.performance_stats['total_traces'] += n_lines
        
        return all_results
        
    def trace_vector_line(self, start_position_re, time_hours, **kwargs):
        """Wrapper bridging single line inputs to the vectorized batch processor."""
        return self.trace_multiple_lines([start_position_re], time_hours, **kwargs)[0]

    def get_performance_stats(self):
        return self.performance_stats

    def benchmark(self, start_points_re, time_hours, n_repeats=3, **kwargs):
        print(f"Starting benchmark with {len(start_points_re)} start points, {n_repeats} repeats each")
        results = {}
        
        for solver in ['euler', 'rk4']:
            config_name = f"solver={solver}"
            print(f"\nTesting configuration: {config_name}")
            
            times = []
            for i in range(n_repeats):
                start_time = time.time()
                traces = self.trace_multiple_lines(start_points_re, time_hours, solver=solver, **kwargs)
                execution_time = time.time() - start_time
                times.append(execution_time)
                print(f"  Run {i+1}/{n_repeats}: {execution_time:.4f} seconds")
            
            avg_time = sum(times) / len(times)
            results[config_name] = {
                "times": times,
                "average": avg_time,
                "throughput": len(start_points_re) / avg_time
            }
            print(f"  Average: {avg_time:.4f} seconds ({results[config_name]['throughput']:.2f} lines/sec)")
        
        return results

def plot_vector_field_traces(traces, tracer, title_suffix="", backend='matplotlib', show_earth=True, interactive=True):
    if backend.lower() == 'plotly':
        return _plot_vector_field_traces_plotly(traces, tracer, title_suffix, show_earth, interactive)
    else:
        return _plot_vector_field_traces_matplotlib(traces, tracer, title_suffix, show_earth)

def _plot_vector_field_traces_matplotlib(traces, tracer, title_suffix="", show_earth=True):
    fig = plt.figure(figsize=(16, 12))
    ax1 = fig.add_subplot(221, projection='3d')
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    
    cmap = get_colormap({'magnetic': 'viridis', 'velocity': 'plasma', 'current': 'inferno', 'electric': 'magma'}.get(tracer.field_type, 'viridis'))
    field_unit = {'magnetic': 'nT', 'velocity': 'km/s', 'current': 'A/m²', 'electric': 'mV/m'}.get(tracer.field_type, 'unit')
    
    if show_earth:
        r_earth = 1.0
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = r_earth * np.outer(np.cos(u), np.sin(v))
        y = r_earth * np.outer(np.sin(u), np.sin(v))
        z = r_earth * np.outer(np.ones(np.size(u)), np.cos(v))
        ax1.plot_surface(x, y, z, color='lightblue', alpha=0.3)
        ax2.add_patch(plt.Circle((0, 0), r_earth, color='lightblue', alpha=0.3))
        ax3.add_patch(plt.Circle((0, 0), r_earth, color='lightblue', alpha=0.3))
    
    all_field_mags = []
    for trace in traces:
        if 'combined' in trace:
            all_field_mags.extend(trace['combined']['field_magnitude'])
        elif 'forward' in trace:
            all_field_mags.extend(trace['forward']['field_magnitude'])
            
    norm = plt.Normalize(vmin=np.min(all_field_mags), vmax=np.max(all_field_mags)) if all_field_mags else plt.Normalize(vmin=0, vmax=1)
    
    for i, trace in enumerate(traces):
        data = trace.get('combined', trace.get('forward'))
        if not data: continue
        
        coords = data['coordinates']
        mags = data['field_magnitude']
        colors = cmap(norm(mags))
        
        ax1.plot(coords[:, 0], coords[:, 1], coords[:, 2], lw=1.5, alpha=0.7)
        for j in range(len(coords)-1):
            ax2.plot(coords[j:j+2, 0], coords[j:j+2, 1], c=colors[j], lw=1.5, alpha=0.7)
            ax3.plot(coords[j:j+2, 0], coords[j:j+2, 2], c=colors[j], lw=1.5, alpha=0.7)
            
        arc_length = np.zeros(len(coords))
        for j in range(1, len(coords)):
            arc_length[j] = arc_length[j-1] + np.linalg.norm(coords[j] - coords[j-1])
        ax4.plot(arc_length, mags, lw=1.5, label=f"Line {i+1}" if len(traces) <= 10 else "")
    
    field_label = tracer.field_type.capitalize()
    ax1.set_title(f"3D {field_label} Field Lines {title_suffix}")
    ax2.set_title(f"XY Plane {title_suffix}")
    ax3.set_title(f"XZ Plane {title_suffix}")
    ax4.set_title(f"{field_label} Field Magnitude Along Lines")
    
    ax1.set_xlabel('X (RE)'); ax1.set_ylabel('Y (RE)'); ax1.set_zlabel('Z (RE)')
    ax2.set_xlabel('X (RE)'); ax2.set_ylabel('Y (RE)')
    ax3.set_xlabel('X (RE)'); ax3.set_ylabel('Z (RE)')
    ax4.set_xlabel('Arc Length (RE)'); ax4.set_ylabel(f'Magnitude ({field_unit})')
    
    ax2.set_aspect('equal')
    ax3.set_aspect('equal')
    if len(traces) <= 10: ax4.legend()
    if all_field_mags:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=[ax1, ax2, ax3], orientation='vertical', pad=0.01, label=f'{field_label} Magnitude ({field_unit})')
        
    plt.tight_layout()
    return fig

def _plot_vector_field_traces_plotly(traces, tracer, title_suffix="", show_earth=True, interactive=True):
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        import plotly.colors as pc
    except ImportError:
        print("Plotly is not installed.")
        return None
        
    fig = make_subplots(rows=2, cols=2, specs=[[{'type': 'scene'}, {'type': 'xy'}], [{'type': 'xy'}, {'type': 'xy'}]])
    field_unit = {'magnetic': 'nT', 'velocity': 'km/s', 'current': 'A/m²', 'electric': 'mV/m'}.get(tracer.field_type, 'unit')
    colorscale = {'magnetic': 'Viridis', 'velocity': 'Plasma', 'current': 'Inferno', 'electric': 'Magma'}.get(tracer.field_type, 'Viridis')
    
    if show_earth:
        theta, phi = np.linspace(0, 2*np.pi, 50), np.linspace(0, np.pi, 50)
        x = np.outer(np.cos(theta), np.sin(phi))
        y = np.outer(np.sin(theta), np.sin(phi))
        z = np.outer(np.ones(np.size(theta)), np.cos(phi))
        fig.add_trace(go.Surface(x=x, y=y, z=z, colorscale=[[0, 'lightblue'], [1, 'lightblue']], opacity=0.3, showscale=False), row=1, col=1)
        
        tc = np.linspace(0, 2*np.pi, 100)
        fig.add_trace(go.Scatter(x=np.cos(tc), y=np.sin(tc), mode='lines', line=dict(color='lightblue', width=2), fill='toself', fillcolor='rgba(173, 216, 230, 0.3)', showlegend=False), row=1, col=2)
        fig.add_trace(go.Scatter(x=np.cos(tc), y=np.zeros(100), mode='lines', line=dict(color='lightblue', width=2), fill='toself', fillcolor='rgba(173, 216, 230, 0.3)', showlegend=False), row=2, col=1)
        
    all_mags = [m for t in traces for m in t.get('combined', t.get('forward', {})).get('field_magnitude', [])]
    f_min, f_max = (np.min(all_mags), np.max(all_mags)) if all_mags else (0, 1)
    
    for i, trace in enumerate(traces):
        data = trace.get('combined', trace.get('forward'))
        if not data: continue
        coords, mags = data['coordinates'], data['field_magnitude']
        
        fig.add_trace(go.Scatter3d(x=coords[:,0], y=coords[:,1], z=coords[:,2], mode='lines', line=dict(color=mags, colorscale=colorscale, width=4, cmin=f_min, cmax=f_max), showlegend=False), row=1, col=1)
        
        line_col = pc.sample_colorscale(colorscale, [(np.mean(mags) - f_min) / (f_max - f_min) if f_max > f_min else 0.5])[0]
        fig.add_trace(go.Scatter(x=coords[:,0], y=coords[:,1], mode='lines', line=dict(color=line_col, width=3), showlegend=False), row=1, col=2)
        fig.add_trace(go.Scatter(x=coords[:,0], y=coords[:,2], mode='lines', line=dict(color=line_col, width=3), showlegend=False), row=2, col=1)
        
        arc = np.zeros(len(coords))
        for j in range(1, len(coords)): arc[j] = arc[j-1] + np.linalg.norm(coords[j] - coords[j-1])
        fig.add_trace(go.Scatter(x=arc, y=mags, mode='lines', line=dict(color=line_col, width=2), showlegend=(len(traces) <= 10)), row=2, col=2)
        
    fig.update_layout(title=f'{tracer.field_type.capitalize()} Field Line Tracing {title_suffix}', height=900)
    fig.update_scenes(xaxis_title='X (RE)', yaxis_title='Y (RE)', zaxis_title='Z (RE)', aspectmode='data')
    return fig

def create_magnetic_tracer(kamodo_object, component_names=None, verbose=False, **kwargs):
    return KamodoVectorFieldTracer(kamodo_object, field_type='magnetic', component_names=component_names or ['B_x', 'B_y', 'B_z'], verbose=verbose, **kwargs)

def create_velocity_tracer(kamodo_object, component_names=None, verbose=False, **kwargs):
    return KamodoVectorFieldTracer(kamodo_object, field_type='velocity', component_names=component_names or ['v_x', 'v_y', 'v_z'], verbose=verbose, **kwargs)

def create_current_tracer(kamodo_object, component_names=None, verbose=False, **kwargs):
    return KamodoVectorFieldTracer(kamodo_object, field_type='current', component_names=component_names or ['J_x', 'J_y', 'J_z'], verbose=verbose, **kwargs)

def create_electric_tracer(kamodo_object, component_names=None, verbose=False, **kwargs):
    return KamodoVectorFieldTracer(kamodo_object, field_type='electric', component_names=component_names or ['E_x', 'E_y', 'E_z'], verbose=verbose, **kwargs)

if __name__ == "__main__":
    print("Optimized Kamodo Vector Field Tracer")

#
### OLD-NEW BREAK
#
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


def benchmark_vectorfield_tracer(kamodo_object, n_lines=10, field_type='magnetic'):
    """
    Comprehensive benchmark of vector field tracer with different optimization settings.
    
    Parameters:
    -----------
    kamodo_object : object
        Kamodo object with field components
    n_lines : int
        Number of field lines to trace for benchmark
    field_type : str
        Type of field to benchmark ('magnetic', 'velocity', etc.)
        
    Returns:
    --------
    results : dict
        Dictionary of benchmark results
    """
    print(f"Starting comprehensive benchmark with {n_lines} field lines")
    
    # Generate random start points in magnetosphere
    rng = np.random.default_rng(42)  # For reproducibility
    theta = rng.uniform(0, 2*np.pi, n_lines)
    phi = rng.uniform(0, np.pi, n_lines)
    r = rng.uniform(3, 10, n_lines)
    
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
    
    # Store results for this field type
    results[field_type] = field_results

    # Print overall summary
    print(f"\n===== BENCHMARK SUMMARY FOR {field_type.upper()} FIELD =====")
    
    baseline = field_results.get('solver=euler', None)
    if baseline is None:
        print("  Baseline configuration not found")
        return results
    
    baseline_time = baseline['execution_time']
    
    # Compare each configuration
    for config_name, config_results in field_results.items():
        speedup = baseline_time / config_results['execution_time'] if config_results['execution_time'] > 0 else 0
        print(f"  {config_name}: {config_results['execution_time']:.4f} sec, "
              f"{speedup:.2f}x speedup relative to Euler, "
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
    print("   # Automatically detects B_x, B_y, B_z components\n")
    
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
    print(f"  Progress bars (tqdm):  {'✓' if TQDM_AVAILABLE else '✗'}")
    
    if not TQDM_AVAILABLE:
        print("\n⚠️  For progress bars, install tqdm:")
        print("    pip install tqdm")
    
    # Show example usage
    print("\n" + "="*50)
    example_usage()
    
    print("\n" + "="*50)
    print("For more examples and detailed documentation,")
    print("see the docstrings and comments in this file.")

