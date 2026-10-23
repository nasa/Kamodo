"""
Kamodo Last Closed Magnetic Field Lines

Functions for finding the last closed magnetic field lines and magnetopause boundary.
All coordinates are in Earth radii (RE).

Updated to use vectorized batch processing:
- Simultaneous bisection for all longitudes
- Explicit RK4 solver enforcing geometric topology
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
    if not trace_result:
        return False
        
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


def _check_closed_batch(tracer, start_pts, time_hours, max_steps=2000):
    """
    Determines if a batch of points are closed field lines using short-circuiting.
    Only traces backward if the forward trace successfully hits Earth.
    """
    if len(start_pts) == 0:
        return np.array([], dtype=bool)
        
    # Trace Forward
    fwd_res = tracer.trace_multiple_lines(
        start_pts, time_hours, direction='forward', solver='rk4', max_steps=max_steps)
    
    # Check which lines successfully hit Earth in the forward direction
    fwd_closed = np.array([('min_altitude' in r['forward']['stop_reason']) for r in fwd_res])
    
    final_closed = np.zeros(len(start_pts), dtype=bool)
    
    # Short-Circuit: Only trace Backward for lines that survived the forward check
    if np.any(fwd_closed):
        bwd_pts = start_pts[fwd_closed]
        bwd_res = tracer.trace_multiple_lines(
            bwd_pts, time_hours, direction='backward', solver='rk4', max_steps=max_steps)
        bwd_closed = np.array([('min_altitude' in r['backward']['stop_reason']) for r in bwd_res])
        
        # Map the surviving backward checks back to the original array
        final_closed[fwd_closed] = bwd_closed
        
    return final_closed


def find_last_closed_field_lines(kamodo_object, time_hours,
    n_traces=36, r_start=3.0, r_max=100.0, tolerance=0.1, show_progress=True, verbose=False):
    """
    Find last closed magnetic field lines and approximate magnetopause boundary.
    
    Uses true vectorized batch processing for simultaneous bisection across all longitudes.
    
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
    tolerance : float, optional (default=0.1)
        Tolerance for bisection method in RE
    show_progress : bool, optional (default=True)
        Whether to show progress bars
    verbose : bool, optional (default=False)
        Whether to print detailed progress information

    Returns:
    --------
    results : dict
        Dictionary containing boundary metrics and performance stats
    """
    if verbose:
        print(f"Finding last closed field lines at t={time_hours:.2f} hours")
        print(f"Using vectorized bisection method with tolerance {tolerance:.3f} RE")
        print(f"Searching from {r_start} to {r_max} RE")
        print(f"Using {n_traces} traces at {360/n_traces:.1f}° intervals")

    start_time = time.time()

    mag_tracer = vft.create_magnetic_tracer(
        kamodo_object,
        verbose=False
    )

    longitudes = np.linspace(0, 360, n_traces, endpoint=False)
    longitude_rad = np.deg2rad(longitudes)
    
    if verbose:
        print(f"Pre-screening inner boundary ({r_start} RE) and outer boundary ({r_max} RE)...")

    # Test Inner Boundary (r_start)
    start_pts = np.column_stack((r_start * np.cos(longitude_rad), r_start * np.sin(longitude_rad), np.zeros(n_traces)))
    start_closed = _check_closed_batch(mag_tracer, start_pts, time_hours, max_steps=2000)

    # Test Outer Boundary (r_max)
    max_pts = np.column_stack((r_max * np.cos(longitude_rad), r_max * np.sin(longitude_rad), np.zeros(n_traces)))
    max_closed = _check_closed_batch(mag_tracer, max_pts, time_hours, max_steps=2000)
    
    # State tracking arrays
    r_low = np.full(n_traces, r_start)
    r_high = np.full(n_traces, r_max)
    
    last_closed_r = np.full(n_traces, np.nan)
    boundary_accuracy = np.full(n_traces, np.nan)
    
    valid_mask = start_closed.copy()
    
    # If a line is already closed at r_max, it escapes our search bounds
    last_closed_r[max_closed] = r_max
    boundary_accuracy[max_closed] = tolerance
    
    # We only bisect lines that start closed but eventually open
    bisect_mask = valid_mask & ~max_closed
    
    max_iterations = int(np.ceil(np.log2((r_max - r_start) / tolerance)))
    
    if verbose:
        print(f"Beginning vectorized bisection (Max iterations: {max_iterations})...")
        
    iter_range = range(max_iterations)
    if show_progress and vft.TQDM_AVAILABLE:
        from tqdm import tqdm
        iter_range = tqdm(iter_range, desc="Bisection search", total=max_iterations)

    for iteration in iter_range:
        if not np.any(bisect_mask):
            break
            
        active_longs = longitude_rad[bisect_mask]
        r_mid = (r_low[bisect_mask] + r_high[bisect_mask]) / 2.0
        
        mid_pts = np.column_stack((r_mid * np.cos(active_longs), r_mid * np.sin(active_longs), np.zeros(len(active_longs))))
        
        # Use short-circuiting in the bisection loop
        mid_closed = _check_closed_batch(mag_tracer, mid_pts, time_hours, max_steps=2000)
        
        active_indices = np.where(bisect_mask)[0]
        closed_idx = active_indices[mid_closed]
        open_idx = active_indices[~mid_closed]
        
        r_low[closed_idx] = r_mid[mid_closed]
        r_high[open_idx] = r_mid[~mid_closed]
        
        bisect_mask = (r_high - r_low) >= tolerance
        bisect_mask = bisect_mask & valid_mask & ~max_closed 

    converged_mask = valid_mask & ~max_closed
    last_closed_r[converged_mask] = (r_low[converged_mask] + r_high[converged_mask]) / 2.0
    boundary_accuracy[converged_mask] = (r_high[converged_mask] - r_low[converged_mask]) / 2.0

    magnetopause_points = []
    for i, r in enumerate(last_closed_r):
        if not np.isnan(r):
            magnetopause_points.append([r * np.cos(longitude_rad[i]), r * np.sin(longitude_rad[i]), 0.0])
    magnetopause_points = np.array(magnetopause_points) if magnetopause_points else np.array([]).reshape(0, 3)

    if verbose:
        print("\nTracing final example field lines for visualization...")

    n_example_traces = min(18, n_traces)
    example_indices = np.linspace(0, n_traces-1, n_example_traces, dtype=int)
    
    closed_start_points = []
    open_start_points = []
    
    for idx in example_indices:
        l_rad = longitude_rad[idx]
        last_r = last_closed_r[idx]
        
        if not np.isnan(last_r):
            r_closed = max(r_start, last_r - 0.5)
            closed_start_points.append([r_closed * np.cos(l_rad), r_closed * np.sin(l_rad), 0.0])
            
            r_open = last_r + 0.5
            open_start_points.append([r_open * np.cos(l_rad), r_open * np.sin(l_rad), 0.0])

    closed_traces = []
    if closed_start_points:
        closed_traces = mag_tracer.trace_multiple_lines(closed_start_points, time_hours, solver='rk4')

    open_traces = []
    if open_start_points:
        open_traces = mag_tracer.trace_multiple_lines(open_start_points, time_hours, solver='rk4')

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
    
    ax1 = fig.add_subplot(221, projection='3d')
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    closed_traces = results['closed_traces']
    open_traces = results['open_traces']
    magnetopause_points = results['magnetopause_points']
    longitudes = results['longitudes']
    last_closed_r = results['last_closed_r']
    summary = results['summary']
    
    if show_earth:
        r_earth = 1.0
        u = np.linspace(0, 2 * np.pi, 50)
        v = np.linspace(0, np.pi, 50)
        x = r_earth * np.outer(np.cos(u), np.sin(v))
        y = r_earth * np.outer(np.sin(u), np.sin(v))
        z = r_earth * np.outer(np.ones(np.size(u)), np.cos(v))
        ax1.plot_surface(x, y, z, color='lightblue', alpha=0.3)
        circle = plt.Circle((0, 0), r_earth, color='lightblue', alpha=0.3)
        ax2.add_patch(circle)

    for trace in closed_traces:
        if 'combined' in trace:
            coords = trace['combined']['coordinates']
        elif 'forward' in trace:
            coords = trace['forward']['coordinates']
        else:
            continue
        ax1.plot(coords[:, 0], coords[:, 1], coords[:, 2], 'b-', alpha=0.7, linewidth=1)
        ax2.plot(coords[:, 0], coords[:, 1], 'b-', alpha=0.7, linewidth=1)

    for trace in open_traces:
        if 'combined' in trace:
            coords = trace['combined']['coordinates']
        elif 'forward' in trace:
            coords = trace['forward']['coordinates']
        else:
            continue
        ax1.plot(coords[:, 0], coords[:, 1], coords[:, 2], 'r-', alpha=0.7, linewidth=1)
        ax2.plot(coords[:, 0], coords[:, 1], 'r-', alpha=0.7, linewidth=1)

    if len(magnetopause_points) > 0:
        boundary_closed = np.vstack([magnetopause_points, magnetopause_points[0]])
        ax1.plot(boundary_closed[:, 0], boundary_closed[:, 1], boundary_closed[:, 2], 
                'g-', linewidth=3, label='Magnetopause')
        ax2.plot(boundary_closed[:, 0], boundary_closed[:, 1], 
                'g-', linewidth=3, label='Magnetopause')
        ax2.scatter(magnetopause_points[:, 0], magnetopause_points[:, 1], 
                   c='green', s=30, alpha=0.7)

    valid_mask = ~np.isnan(last_closed_r)
    if np.any(valid_mask):
        ax3.plot(longitudes[valid_mask], last_closed_r[valid_mask], 'go-', 
                linewidth=2, markersize=4)
        ax3.fill_between(longitudes[valid_mask], 
                        last_closed_r[valid_mask] - results['boundary_accuracy'][valid_mask],
                        last_closed_r[valid_mask] + results['boundary_accuracy'][valid_mask],
                        alpha=0.3, color='green')

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

    ax1.set_xlabel('X (RE)'); ax1.set_ylabel('Y (RE)'); ax1.set_zlabel('Z (RE)')
    ax1.set_title('3D View: Last Closed Magnetic Field Lines')

    ax2.set_xlabel('X (RE)'); ax2.set_ylabel('Y (RE)')
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
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{'type': 'scene'}, {'type': 'xy'}],
               [{'type': 'xy'}, {'type': 'xy'}]],
        subplot_titles=('3D View: Last Closed Magnetic Field Lines',
                       'XY Plane View',
                       'Magnetopause Boundary Profile',
                       'Performance Summary')
    )

    closed_traces = results['closed_traces']
    open_traces = results['open_traces']
    magnetopause_points = results['magnetopause_points']
    longitudes = results['longitudes']
    last_closed_r = results['last_closed_r']
    summary = results['summary']

    if show_earth:
        r_earth = 1.0
        theta = np.linspace(0, 2 * np.pi, 50)
        phi = np.linspace(0, np.pi, 50)
        x = r_earth * np.outer(np.cos(theta), np.sin(phi))
        y = r_earth * np.outer(np.sin(theta), np.sin(phi))
        z = r_earth * np.outer(np.ones(np.size(theta)), np.cos(phi))
        
        fig.add_trace(go.Surface(x=x, y=y, z=z, colorscale=[[0, 'lightblue'], [1, 'lightblue']], opacity=0.3, showscale=False, name='Earth'), row=1, col=1)
        
        theta_circle = np.linspace(0, 2*np.pi, 100)
        x_circle = r_earth * np.cos(theta_circle)
        y_circle = r_earth * np.sin(theta_circle)
        
        fig.add_trace(go.Scatter(x=x_circle, y=y_circle, mode='lines', line=dict(color='lightblue', width=2), fill='toself', fillcolor='rgba(173, 216, 230, 0.3)', showlegend=False, name='Earth'), row=1, col=2)

    for i, trace in enumerate(closed_traces):
        if 'combined' in trace: coords = trace['combined']['coordinates']
        elif 'forward' in trace: coords = trace['forward']['coordinates']
        else: continue
        
        fig.add_trace(go.Scatter3d(x=coords[:, 0], y=coords[:, 1], z=coords[:, 2], mode='lines', line=dict(color='blue', width=2), showlegend=(i==0), name='Closed field lines'), row=1, col=1)
        fig.add_trace(go.Scatter(x=coords[:, 0], y=coords[:, 1], mode='lines', line=dict(color='blue', width=2), showlegend=False), row=1, col=2)

    for i, trace in enumerate(open_traces):
        if 'combined' in trace: coords = trace['combined']['coordinates']
        elif 'forward' in trace: coords = trace['forward']['coordinates']
        else: continue
        
        fig.add_trace(go.Scatter3d(x=coords[:, 0], y=coords[:, 1], z=coords[:, 2], mode='lines', line=dict(color='red', width=2), showlegend=(i==0), name='Open field lines'), row=1, col=1)
        fig.add_trace(go.Scatter(x=coords[:, 0], y=coords[:, 1], mode='lines', line=dict(color='red', width=2), showlegend=False), row=1, col=2)

    if len(magnetopause_points) > 0:
        boundary_closed = np.vstack([magnetopause_points, magnetopause_points[0]])
        fig.add_trace(go.Scatter3d(x=boundary_closed[:, 0], y=boundary_closed[:, 1], z=boundary_closed[:, 2], mode='lines', line=dict(color='green', width=4), name='Magnetopause'), row=1, col=1)
        fig.add_trace(go.Scatter(x=boundary_closed[:, 0], y=boundary_closed[:, 1], mode='lines+markers', line=dict(color='green', width=4), marker=dict(color='green', size=6), showlegend=False), row=1, col=2)

    valid_mask = ~np.isnan(last_closed_r)
    if np.any(valid_mask):
        fig.add_trace(go.Scatter(x=longitudes[valid_mask], y=last_closed_r[valid_mask], mode='lines+markers', line=dict(color='green', width=3), marker=dict(color='green', size=6), name='Boundary radius', showlegend=False), row=2, col=1)
        
        upper_bound = last_closed_r[valid_mask] + results['boundary_accuracy'][valid_mask]
        lower_bound = last_closed_r[valid_mask] - results['boundary_accuracy'][valid_mask]
        
        fig.add_trace(go.Scatter(x=np.concatenate([longitudes[valid_mask], longitudes[valid_mask][::-1]]), y=np.concatenate([upper_bound, lower_bound[::-1]]), fill='toself', fillcolor='rgba(0,255,0,0.3)', line=dict(color='rgba(255,255,255,0)'), showlegend=False), row=2, col=1)

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
    
    fig.add_annotation(text=summary_text, xref="paper", yref="paper", xshift=0, yshift=100, x=0.5, y=1.0, xanchor='center', yanchor='top', showarrow=False, font=dict(size=12, family="monospace"), row=2, col=2)

    fig.update_scenes(xaxis_title='X (RE)', yaxis_title='Y (RE)', zaxis_title='Z (RE)', aspectmode='data', row=1, col=1)
    fig.update_xaxes(title_text='X (RE)', scaleanchor='y', scaleratio=1, row=1, col=2)
    fig.update_yaxes(title_text='Y (RE)', row=1, col=2)
    fig.update_xaxes(title_text='Longitude (degrees)', range=[0, 360], row=2, col=1)
    fig.update_yaxes(title_text='Last Closed Radius (RE)', row=2, col=1)
    fig.update_xaxes(showticklabels=False, showgrid=False, zeroline=False, row=2, col=2)
    fig.update_yaxes(showticklabels=False, showgrid=False, zeroline=False, row=2, col=2)
    
    fig.update_layout(title='Last Closed Magnetic Field Lines and Magnetopause Boundary', height=800, width=1200, hovermode='closest' if interactive else False, legend=dict(x=0.01, y=0.99, bgcolor='rgba(255,255,255,0.7)', bordercolor='rgba(0,0,0,0.5)', borderwidth=1))
    
    return fig


def find_last_closed_field_lines_auto_bounds(ko, time_hours, n_traces=36,
    r_start=2.0, r_max_factor=0.8, tolerance=0.01, show_progress=True, verbose=False):
    """
    Find last closed magnetic field lines with automatic bound detection.
    """
    test_tracer = vft.create_magnetic_tracer(
        ko,
        verbose=False
    )

    try:
        bounds = {'x_range': [-30, 30], 'y_range': [-30, 30], 'z_range': [-30, 30]}
        test_radius = 25.0
        test_points = []
        
        for r in [5.0, 10.0, 15.0, 20.0, test_radius]:
            for theta in [0, np.pi/2, np.pi, 3*np.pi/2]:
                test_points.append([r * np.cos(theta), r * np.sin(theta), 0.0])
        
        valid_points = []
        for point in test_points:
            try:
                test_trace = test_tracer.trace_vector_line(point, time_hours, solver='rk4')
                if test_trace is not None:
                    valid_points.append(point)
            except Exception:
                pass 
                
        if valid_points:
            valid_points = np.array(valid_points)
            max_r = np.max(np.sqrt(np.sum(valid_points[:, :2]**2, axis=1)))
            bounds = {
                'x_range': [-max_r, max_r],
                'y_range': [-max_r, max_r],
                'z_range': [-10, 10] 
            }
            if verbose:
                print(f"Detected model bounds: r_max={max_r:.1f} RE")
        else:
            if verbose:
                print(f"Could not detect model bounds, using defaults")
    except Exception as e:
        if verbose:
            print(f"Error detecting bounds: {e}, using defaults")
    
    x_extent = bounds['x_range'][1] - bounds['x_range'][0]
    y_extent = bounds['y_range'][1] - bounds['y_range'][0]
    model_radius = min(x_extent, y_extent) / 2
    r_max = min(r_max_factor * model_radius, 20.0)
    
    if verbose:
        print(f"Using r_max = {r_max:.1f} RE (model_radius = {model_radius:.1f} RE)")
    
    results = find_last_closed_field_lines(
        ko, time_hours,
        n_traces=n_traces,
        r_start=r_start,
        r_max=r_max,
        tolerance=tolerance,
        show_progress=show_progress,
        verbose=verbose
    )
    
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
    """
    magnetopause_points = results['magnetopause_points']
    longitudes = results['longitudes']
    last_closed_r = results['last_closed_r']
    
    if len(magnetopause_points) == 0 or not np.any(~np.isnan(last_closed_r)):
        return {
            'status': 'error',
            'message': 'No valid magnetopause boundary points found'
        }
    
    valid_mask = ~np.isnan(last_closed_r)
    valid_longitudes = longitudes[valid_mask]
    valid_radii = last_closed_r[valid_mask]
    valid_points = magnetopause_points
    
    mean_radius = np.mean(valid_radii)
    std_radius = np.std(valid_radii)
    min_radius = np.min(valid_radii)
    min_long = valid_longitudes[np.argmin(valid_radii)]
    max_radius = np.max(valid_radii)
    max_long = valid_longitudes[np.argmax(valid_radii)]
    
    day_mask = (valid_longitudes >= 270) | (valid_longitudes <= 90)
    night_mask = ~day_mask
    
    day_mean = np.mean(valid_radii[day_mask]) if np.any(day_mask) else np.nan
    night_mean = np.mean(valid_radii[night_mask]) if np.any(night_mask) else np.nan
    day_night_ratio = day_mean / night_mean if not np.isnan(night_mean) and night_mean > 0 else np.nan
    
    dawn_mask = (valid_longitudes >= 0) & (valid_longitudes <= 180)
    dusk_mask = ~dawn_mask
    
    dawn_mean = np.mean(valid_radii[dawn_mask]) if np.any(dawn_mask) else np.nan
    dusk_mean = np.mean(valid_radii[dusk_mask]) if np.any(dusk_mask) else np.nan
    dawn_dusk_ratio = dawn_mean / dusk_mean if not np.isnan(dusk_mean) and dusk_mean > 0 else np.nan
    
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
                               tolerance_options=[0.2, 0.05],
                               verbose=True):
    """
    Benchmark different tolerance configurations of the last closed field line finder.
    """
    print(f"Starting benchmark with {n_traces} longitude slices")
    
    results = {}
    
    for tolerance in tolerance_options:
        config_name = f"tolerance={tolerance}"
        print(f"\nTesting configuration: {config_name}")
        
        start_time = time.time()
        
        result = find_last_closed_field_lines(
            ko, time_hours,
            n_traces=n_traces,
            r_start=2.0,
            r_max=12.0,
            tolerance=tolerance,
            verbose=verbose
        )
        
        execution_time = time.time() - start_time
        
        results[config_name] = {
            'execution_time': execution_time,
            'valid_boundaries': result['summary']['n_valid_boundaries'],
            'mean_radius': result['summary']['mean_boundary_radius'],
            'longitudes_per_second': n_traces / execution_time
        }
        
        print(f"Execution time: {execution_time:.2f} seconds")
        print(f"Performance: {n_traces / execution_time:.2f} longitudes/second")
        print(f"Valid boundaries: {result['summary']['n_valid_boundaries']}/{n_traces}")
    
    if results:
        print("\n===== BENCHMARK SUMMARY =====")
        baseline = results.get(f"tolerance={tolerance_options[0]}", None)
        baseline_time = baseline['execution_time'] if baseline else None
        
        for config_name, data in results.items():
            speedup = baseline_time / data['execution_time'] if baseline_time else 1.0
            print(f"{config_name}: {data['execution_time']:.2f} sec, "
                  f"{speedup:.2f}x relative speed, "
                  f"{data['longitudes_per_second']:.2f} longitudes/sec")
    
    return results


def example_usage(ko, time_hours=12.0):
    """
    Example of finding and plotting last closed field lines.
    """
    print("Example: Finding Last Closed Field Lines")
    print("=" * 40)
    
    results = find_last_closed_field_lines_auto_bounds(
        ko, time_hours,
        n_traces=36,
        r_start=2.0,
        verbose=True,
        show_progress=True
    )
    
    shape_analysis = analyze_magnetopause_shape(results)
    
    print("\nCreating plots...")
    fig_mpl = plot_last_closed_field_lines(results, backend='matplotlib')
    plt.show()
    
    if PLOTLY_AVAILABLE:
        fig_plotly = plot_last_closed_field_lines(results, backend='plotly')
        try:
            from IPython import get_ipython
            if get_ipython() is not None:
                fig_plotly.show()
        except ImportError:
            print("To view interactive plot, call: fig_plotly.show()")
    
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
    print("Kamodo Last Closed Field Line Finder")
    print("To use this module, import it and call its functions.")
    print("For an example, use:")
    print("  results = find_last_closed_field_lines(kamodo_object, time_hours)")
    print("  fig = plot_last_closed_field_lines(results)")
    
    print("\nAvailable tools:")
    print(f"  Progress bars (tqdm):  {'✓' if vft.TQDM_AVAILABLE else '✗'}")
