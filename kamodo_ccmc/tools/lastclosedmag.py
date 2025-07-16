"""
Kamodo Last Closed Magnetic Field Lines

Functions for finding the last closed magnetic field lines and magnetopause boundary.
All coordinates are in Earth radii (RE).
"""

import numpy as np
import matplotlib.pyplot as plt
try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

# Import the vector field tracer
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

    # Count how many ends hit Earth's surface
    earth_hits = sum(1 for reason in stop_reasons if 'Earth_surface' in reason)

    # Field line is closed if both ends hit Earth
    is_closed = earth_hits >= 2

    return is_closed


def _find_boundary_bisection(mag_tracer, lon_rad, time_hours, r_start, r_max,
                            tolerance, max_iter, initial_step, max_steps, step_size_re,
                            verbose=False):
    """
    Use bisection method to find the last closed field line boundary accurately.

    Parameters:
    -----------
    mag_tracer : KamodoVectorFieldTracer
        Magnetic field tracer
    lon_rad : float
        Longitude in radians
    time_hours : float
        Time in hours
    r_start, r_max : float
        Search range in RE
    tolerance : float
        Desired accuracy in RE
    max_iter : int
        Maximum bisection iterations
    initial_step : float
        Initial step for coarse search
    max_steps : int
        Max integration steps per trace
    step_size_re : float
        Integration step size
    verbose : bool
        Print detailed progress

    Returns:
    --------
    boundary_r : float
        Last closed field line radius (NaN if not found)
    final_accuracy : float
        Achieved accuracy
    """
    # Step 1: Coarse search to bracket the boundary
    if verbose:
        print(f"      Coarse search from {r_start:.1f} to {r_max:.1f} RE...")

    r_closed = None  # Last known closed radius
    r_open = None    # First known open radius

    # Search with initial step size to bracket the boundary
    r_test = r_start
    while r_test <= r_max:
        x_test = r_test * np.cos(lon_rad)
        y_test = r_test * np.sin(lon_rad)
        z_test = 0.0

        try:
            trace_result = mag_tracer.trace_vector_line(
                (x_test, y_test, z_test), time_hours,
                max_steps=max_steps,
                step_size_re=step_size_re,
                direction='both',
                adaptive_step=True
            )

            is_closed = _is_field_line_closed(trace_result)

            if is_closed:
                r_closed = r_test
                if verbose:
                    print(f"        R={r_test:.2f} RE: CLOSED")
            else:
                r_open = r_test
                if verbose:
                    print(f"        R={r_test:.2f} RE: OPEN")
                break  # Found first open field line

        except Exception as e:
            if verbose:
                print(f"        R={r_test:.2f} RE: ERROR ({e})")
            break

        r_test += initial_step

    # Check if we found a valid bracket
    if r_closed is None or r_open is None:
        if verbose:
            print(f"      Could not bracket boundary (closed={r_closed}, open={r_open})")
        return np.nan, np.nan

    # Step 2: Bisection method for accurate boundary finding
    if verbose:
        print(f"      Bisection between {r_closed:.3f} and {r_open:.3f} RE...")

    r_low = r_closed    # Known closed
    r_high = r_open     # Known open

    for iteration in range(max_iter):
        r_mid = (r_low + r_high) / 2.0
        current_accuracy = (r_high - r_low) / 2.0

        if current_accuracy <= tolerance:
            if verbose:
                print(f"      Converged after {iteration} iterations: "
                      f"{r_mid:.4f} ± {current_accuracy:.4f} RE")
            return r_mid, current_accuracy

        # Test the midpoint
        x_mid = r_mid * np.cos(lon_rad)
        y_mid = r_mid * np.sin(lon_rad)
        z_mid = 0.0

        try:
            trace_result = mag_tracer.trace_vector_line(
                (x_mid, y_mid, z_mid), time_hours,
                max_steps=max_steps,
                step_size_re=step_size_re,
                direction='both',
                adaptive_step=True
            )

            is_closed = _is_field_line_closed(trace_result)

            if is_closed:
                r_low = r_mid  # Midpoint is closed, move lower bound up
                if verbose:
                    print(f"        Iter {iteration+1}: R={r_mid:.4f} RE CLOSED, "
                          f"bracket=[{r_low:.4f}, {r_high:.4f}], acc=±{current_accuracy:.4f}")
            else:
                r_high = r_mid  # Midpoint is open, move upper bound down
                if verbose:
                    print(f"        Iter {iteration+1}: R={r_mid:.4f} RE OPEN, "
                          f"bracket=[{r_low:.4f}, {r_high:.4f}], acc=±{current_accuracy:.4f}")

        except Exception as e:
            if verbose:
                print(f"        Iter {iteration+1}: R={r_mid:.4f} RE ERROR: {e}")
            # On error, assume open (conservative approach)
            r_high = r_mid

    # Return best estimate after max iterations
    final_r = (r_low + r_high) / 2.0
    final_accuracy = (r_high - r_low) / 2.0

    if verbose:
        print(f"      Max iterations reached: {final_r:.4f} ± {final_accuracy:.4f} RE")

    return final_r, final_accuracy


def find_last_closed_field_lines(ko, time_hours, n_traces=18, r_start=2.0, r_max=20.0,
                                tolerance=0.02, max_bisection_iter=10, initial_step=4.0,
                                max_steps=5000, step_size_re=0.1,
                                coord_system='GSM', verbose=True):
    """
    Find the last closed magnetic field lines using bisection method for accurate boundary detection.
    Seeds points are placed in the equatorial plane at evenly spaced longitudes.

    Parameters:
    -----------
    ko : Kamodo object
        Kamodo object containing magnetic field data
    time_hours : float
        Time in hours since midnight of first day of data
    n_traces : int
        Number of traces (default 18 for 20-degree spacing)
    r_start : float
        Starting radial distance in RE (default 2.0 RE)
    r_max : float
        Maximum radial distance to search in RE (default 20.0 RE)
    tolerance : float
        Tolerance for bisection method in RE (default 0.02 RE)
    max_bisection_iter : int
        Maximum iterations for bisection method (default 10)
    initial_step : float
        Initial step size for coarse search in RE (default 4.0 RE)
    max_steps : int
        Maximum integration steps per trace
    step_size_re : float
        Integration step size in RE
    coord_system : str
        Coordinate system ('GSM', 'GSE', 'SM', etc.)
    verbose : bool
        Whether to print progress information

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
    """
    if verbose:
        print(f"Finding last closed field lines at t={time_hours:.2f} hours")
        print(f"Using bisection method with tolerance {tolerance:.3f} RE")
        print(f"Searching from {r_start} to {r_max} RE")
        print(f"Using {n_traces} traces at {360/n_traces:.1f}° intervals")

    # Create magnetic field tracer
    mag_tracer = vft.KamodoVectorFieldTracer(
        ko, vector_components=('B_x', 'B_y', 'B_z'),
        coord_system=coord_system, field_type='magnetic')

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

    # For each longitude direction
    for i, lon_deg in enumerate(longitudes):
        if verbose and (i + 1) % progress_interval == 0:
            print(f"  Processing longitude {i+1}/{n_traces} ({lon_deg:.1f}°)")

        lon_rad = np.radians(lon_deg)

        # Find boundary using bisection method
        boundary_r, final_accuracy = _find_boundary_bisection(
            mag_tracer, lon_rad, time_hours, r_start, r_max,
            tolerance, max_bisection_iter, initial_step,
            max_steps, step_size_re, verbose and i < 3  # Only show details for first few
        )

        if not np.isnan(boundary_r):
            last_closed_r[i] = boundary_r
            boundary_accuracy[i] = final_accuracy

            # Get the final closed field line trace
            x_closed = boundary_r * np.cos(lon_rad)
            y_closed = boundary_r * np.sin(lon_rad)
            z_closed = 0.0
            closed_start_point = (x_closed, y_closed, z_closed)

            try:
                closed_trace = mag_tracer.trace_vector_line(
                    closed_start_point, time_hours,
                    max_steps=max_steps,
                    step_size_re=step_size_re,
                    direction='both',
                    adaptive_step=True
                )

                closed_traces.append({
                    'trace': closed_trace,
                    'longitude': lon_deg,
                    'radius': boundary_r,
                    'start_point': closed_start_point,
                    'accuracy': final_accuracy
                })

                # Get the first open field line trace (slightly beyond boundary)
                open_r = boundary_r + 2 * tolerance  # Slightly beyond closed boundary
                x_open = open_r * np.cos(lon_rad)
                y_open = open_r * np.sin(lon_rad)
                z_open = 0.0
                open_start_point = (x_open, y_open, z_open)

                open_trace = mag_tracer.trace_vector_line(
                    open_start_point, time_hours,
                    max_steps=max_steps,
                    step_size_re=step_size_re,
                    direction='both',
                    adaptive_step=True
                )

                open_traces.append({
                    'trace': open_trace,
                    'longitude': lon_deg,
                    'radius': open_r,
                    'start_point': open_start_point
                })

                # Add magnetopause point (use the open field line starting point)
                magnetopause_points.append([x_open, y_open, z_open])

                if verbose:
                    print(f"    Lon {lon_deg:.1f}°: Last closed at {boundary_r:.3f} RE "
                          f"(±{final_accuracy:.3f} RE)")

            except Exception as e:
                if verbose:
                    print(f"    Error getting final traces at lon={lon_deg:.1f}°: {e}")
        else:
            if verbose:
                print(f"    Lon {lon_deg:.1f}°: No boundary found in search range")

    # Create magnetopause boundary estimate
    magnetopause_points = (np.array(magnetopause_points)
                           if magnetopause_points else np.array([]).reshape(0, 3))

    # Summary statistics
    valid_last_closed = last_closed_r[~np.isnan(last_closed_r)]
    valid_accuracy = boundary_accuracy[~np.isnan(boundary_accuracy)]

    if verbose:
        print(f"\nSummary:")
        print(f"  Found closed field lines at {len(valid_last_closed)}/{n_traces} longitudes")
        if len(valid_last_closed) > 0:
            print(f"  Last closed radius: mean={np.mean(valid_last_closed):.3f} RE, "
                  f"range=[{np.min(valid_last_closed):.3f}, {np.max(valid_last_closed):.3f}] RE")
            print(f"  Average accuracy: {np.mean(valid_accuracy):.4f} RE")
            print(f"  Magnetopause points: {len(magnetopause_points)}")
        else:
            print(f"  No closed field lines found in search range")
        print(f"  Number of traces computed: {mag_tracer.report_trace_usage(value=True)}")

    results = {
        'closed_traces': closed_traces,
        'open_traces': open_traces,
        'last_closed_r': last_closed_r,
        'longitudes': longitudes,
        'magnetopause_points': magnetopause_points,
        'boundary_accuracy': boundary_accuracy,
        'summary': {
            'n_closed': len(valid_last_closed),
            'n_total': n_traces,
            'mean_radius': np.mean(valid_last_closed) if len(valid_last_closed) > 0 else np.nan,
            'min_radius': np.min(valid_last_closed) if len(valid_last_closed) > 0 else np.nan,
            'max_radius': np.max(valid_last_closed) if len(valid_last_closed) > 0 else np.nan,
            'mean_accuracy': np.mean(valid_accuracy) if len(valid_accuracy) > 0 else np.nan,
            'time_hours': time_hours,
            'coord_system': coord_system,
            'tolerance': tolerance,
            'trace_calls': mag_tracer.report_trace_usage(value=True)
        }
    }

    return results


def find_last_closed_field_lines_auto_bounds(ko, time_hours, n_traces=36,
                                            r_start=2.0, r_max_factor=0.8,
                                            tolerance=0.01, max_bisection_iter=10,
                                            initial_step=1.0, max_steps=5000,
                                            step_size_re=0.1, coord_system='GSM',
                                            verbose=True):
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
    verbose : bool
        Whether to print progress information
    
    Returns:
    --------
    results : dict
        Dictionary containing trace results and model bounds information
    """
    # Get model bounds automatically
    vector_components = ('B_x', 'B_y', 'B_z')
    bounds = vft.get_model_bounds(ko, vector_components, verbose=verbose)
    
    # Calculate search parameters based on model bounds
    max_extent = bounds['max_extent']
    r_max = min(r_max_factor * max_extent, max_extent * 0.9)  # Don't go to edge
    
    if verbose:
        print(f"\nAutomatic search parameter calculation:")
        print(f"  Model max extent: {max_extent:.1f} RE")
        print(f"  Search range: {r_start:.1f} to {r_max:.1f} RE")
        print(f"  r_start: {r_start:.1f} RE (fixed starting radius)")
        print(f"  r_max_factor: {r_max_factor} -> r_max: {r_max:.1f} RE")
    
    # Use the main function with calculated bounds
    results = find_last_closed_field_lines(
        ko, time_hours, n_traces=n_traces,
        r_start=r_start, r_max=r_max,
        tolerance=tolerance, max_bisection_iter=max_bisection_iter,
        initial_step=initial_step, max_steps=max_steps,
        step_size_re=step_size_re, coord_system=coord_system,
        verbose=verbose
    )
    
    # Add bounds information to results
    results['model_bounds'] = bounds
    results['search_parameters'] = {
        'r_start': r_start,
        'r_max': r_max,
        'r_max_factor': r_max_factor
    }
    
    return results


def plot_last_closed_field_lines(results, backend='matplotlib', show_magnetopause=True):
    """
    Plot the last closed field lines and magnetopause boundary.
    
    Parameters:
    -----------
    results : dict
        Results from find_last_closed_field_lines
    backend : str
        'matplotlib' or 'plotly'
    show_magnetopause : bool
        Whether to show estimated magnetopause boundary
    
    Returns:
    --------
    fig : matplotlib or plotly figure
        The generated plot
    """
    if backend.lower() == 'plotly':
        if not PLOTLY_AVAILABLE:
            print("Plotly is not available. Install with: pip install plotly")
            backend = 'matplotlib'
        else:
            return _plot_last_closed_field_lines_plotly(results, show_magnetopause)
    
    return _plot_last_closed_field_lines_matplotlib(results, show_magnetopause)


def _plot_last_closed_field_lines_matplotlib(results, show_magnetopause=True):
    """Create matplotlib plot of last closed field lines with aspect ratio 1."""
    
    fig = plt.figure(figsize=(18, 12))
    
    # Create subplots
    ax1 = fig.add_subplot(221, projection='3d')  # 3D view
    ax2 = fig.add_subplot(222)  # XY plane (equatorial)
    ax3 = fig.add_subplot(223)  # Radial distance vs longitude
    ax4 = fig.add_subplot(224)  # XZ plane (meridional)
    
    closed_traces = results['closed_traces']
    open_traces = results['open_traces']
    magnetopause_points = results['magnetopause_points']
    longitudes = results['longitudes']
    last_closed_r = results['last_closed_r']
    
    # Collect all coordinates for aspect ratio calculation
    all_coords = []
    
    # Plot closed field lines
    for i, trace_data in enumerate(closed_traces):
        if trace_data is None:
            continue
            
        trace = trace_data['trace']
        trace_coords = trace.get('combined', trace.get('forward', {})).get('coordinates', np.array([]))
        
        if len(trace_coords) < 2:
            continue
        
        all_coords.extend(trace_coords)
        x, y, z = trace_coords.T
        color = plt.cm.viridis(i / len(closed_traces))
        
        # 3D plot
        ax1.plot(x, y, z, color=color, alpha=0.8, linewidth=1.5)
        ax1.scatter(*trace_data['start_point'], color=color, s=50, marker='o')
        
        # XY projection
        ax2.plot(x, y, color=color, alpha=0.8, linewidth=1.5)
        ax2.scatter(*trace_data['start_point'][:2], color=color, s=50, marker='o')
        
        # XZ projection  
        ax4.plot(x, z, color=color, alpha=0.8, linewidth=1.5)
    
    # Plot open field lines (first open field line at each longitude)
    for i, trace_data in enumerate(open_traces):
        if trace_data is None:
            continue
            
        trace = trace_data['trace']
        trace_coords = trace.get('combined', trace.get('forward', {})).get('coordinates', np.array([]))
        
        if len(trace_coords) < 2:
            continue
        
        all_coords.extend(trace_coords)
        x, y, z = trace_coords.T
        
        # Plot in red to distinguish from closed field lines
        ax1.plot(x, y, z, color='red', alpha=0.6, linewidth=1, linestyle='--')
        ax2.plot(x, y, color='red', alpha=0.6, linewidth=1, linestyle='--')
        ax4.plot(x, z, color='red', alpha=0.6, linewidth=1, linestyle='--')
    
    # Plot magnetopause boundary estimate
    if show_magnetopause and len(magnetopause_points) > 0:
        mp_x, mp_y, mp_z = magnetopause_points.T
        
        # 3D scatter
        ax1.scatter(mp_x, mp_y, mp_z, color='orange', s=100, marker='s', 
                   label='Magnetopause', alpha=0.8)
        
        # XY projection with curve fit
        ax2.scatter(mp_x, mp_y, color='orange', s=100, marker='s', alpha=0.8)
        
        # Try to fit a smooth curve through magnetopause points
        if len(mp_x) > 3:
            try:
                # Sort by angle for smooth curve
                angles = np.arctan2(mp_y, mp_x)
                sort_idx = np.argsort(angles)
                
                # Add first point at end to close the curve
                mp_x_sorted = np.append(mp_x[sort_idx], mp_x[sort_idx[0]])
                mp_y_sorted = np.append(mp_y[sort_idx], mp_y[sort_idx[0]])
                
                ax2.plot(mp_x_sorted, mp_y_sorted, 'orange', linewidth=3, 
                        alpha=0.7, label='Magnetopause boundary')
            except Exception:
                pass  # Silently handle any curve fitting errors
    
    # Plot radial distance vs longitude with error bars
    valid_mask = ~np.isnan(last_closed_r)
    if np.any(valid_mask):
        boundary_accuracy = results.get('boundary_accuracy', np.full_like(last_closed_r, 0.01))
        ax3.errorbar(longitudes[valid_mask], last_closed_r[valid_mask], 
                    yerr=boundary_accuracy[valid_mask],
                    fmt='bo-', linewidth=2, markersize=6, capsize=3,
                    label='Last closed radius')
        ax3.set_xlabel('Longitude [degrees]')
        ax3.set_ylabel('Last Closed Radius [RE]')
        ax3.set_title('Last Closed Field Line Radius vs Longitude')
        ax3.grid(True, alpha=0.3)
        ax3.legend()
    
    # Add Earth to 3D and 2D plots
    # 3D Earth
    u = np.linspace(0, 2 * np.pi, 20)
    v = np.linspace(0, np.pi, 20)
    x_earth = np.outer(np.cos(u), np.sin(v))
    y_earth = np.outer(np.sin(u), np.sin(v))
    z_earth = np.outer(np.ones(np.size(u)), np.cos(v))
    ax1.plot_surface(x_earth, y_earth, z_earth, alpha=0.3, color='blue')
    
    # Set aspect ratio to 1 for 3D plot
    if all_coords:
        all_coords = np.array(all_coords)
        max_range = np.array([all_coords[:, 0].max()-all_coords[:, 0].min(),
                             all_coords[:, 1].max()-all_coords[:, 1].min(),
                             all_coords[:, 2].max()-all_coords[:, 2].min()]).max() / 2.0
        mid_x = (all_coords[:, 0].max()+all_coords[:, 0].min()) * 0.5
        mid_y = (all_coords[:, 1].max()+all_coords[:, 1].min()) * 0.5
        mid_z = (all_coords[:, 2].max()+all_coords[:, 2].min()) * 0.5
        ax1.set_xlim(mid_x - max_range, mid_x + max_range)
        ax1.set_ylim(mid_y - max_range, mid_y + max_range)
        ax1.set_zlim(mid_z - max_range, mid_z + max_range)
    
    # 2D Earth circles with aspect ratio 1
    theta = np.linspace(0, 2*np.pi, 100)
    for ax in [ax2, ax4]:
        ax.plot(np.cos(theta), np.sin(theta), 'b-', linewidth=2, alpha=0.7)
        ax.set_aspect('equal')  # Aspect ratio 1
        ax.grid(True, alpha=0.3)
    
    # Set labels and titles
    coord_sys = results['summary']['coord_system']
    time_hours = results['summary']['time_hours']
    tolerance = results['summary'].get('tolerance', 0.01)
    
    ax1.set_xlabel(f'X ({coord_sys}) [RE]')
    ax1.set_ylabel(f'Y ({coord_sys}) [RE]')
    ax1.set_zlabel(f'Z ({coord_sys}) [RE]')
    ax1.set_title(f'Last Closed Field Lines 3D (t={time_hours:.1f}h, tol={tolerance:.3f})')
    
    ax2.set_xlabel(f'X ({coord_sys}) [RE]')
    ax2.set_ylabel(f'Y ({coord_sys}) [RE]')
    ax2.set_title('Equatorial Plane View')
    
    ax4.set_xlabel(f'X ({coord_sys}) [RE]')
    ax4.set_ylabel(f'Z ({coord_sys}) [RE]')
    ax4.set_title('Meridional Plane View')
    
    plt.tight_layout()
    return fig


def _plot_last_closed_field_lines_plotly(results, show_magnetopause=True):
    """Create Plotly plot of last closed field lines with aspect ratio 1."""
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{"type": "scene", "rowspan": 2}, {"type": "xy"}],
               [None, {"type": "xy"}]],
        subplot_titles=[
            'Last Closed Field Lines 3D',
            'Equatorial Plane View',
            'Last Closed Radius vs Longitude'
        ],
        horizontal_spacing=0.1,
        vertical_spacing=0.1
    )
    
    closed_traces = results['closed_traces']
    open_traces = results['open_traces']
    magnetopause_points = results['magnetopause_points']
    longitudes = results['longitudes']
    last_closed_r = results['last_closed_r']
    
    # Plot closed field lines
    for i, trace_data in enumerate(closed_traces):
        if trace_data is None:
            continue
            
        trace = trace_data['trace']
        trace_coords = trace.get('combined', trace.get('forward', {})).get('coordinates', np.array([]))
        
        if len(trace_coords) < 2:
            continue
        
        x, y, z = trace_coords.T
        color = px.colors.sample_colorscale('Viridis', [i / len(closed_traces)])[0]
        accuracy = trace_data.get('accuracy', 0.01)
        
        # 3D plot
        fig.add_trace(
            go.Scatter3d(
                x=x, y=y, z=z,
                mode='lines',
                line=dict(color=color, width=4),
                name=f'Closed {trace_data["longitude"]:.0f}°',
                showlegend=False,
                hovertemplate=(f'Longitude: {trace_data["longitude"]:.1f}°<br>'
                             f'Radius: {trace_data["radius"]:.3f} RE<br>'
                             f'Accuracy: ±{accuracy:.3f} RE<extra></extra>')
            ),
            row=1, col=1
        )
        
        # Starting point
        fig.add_trace(
            go.Scatter3d(
                x=[trace_data['start_point'][0]], 
                y=[trace_data['start_point'][1]], 
                z=[trace_data['start_point'][2]],
                mode='markers',
                marker=dict(color=color, size=8),
                name=f'Start {trace_data["longitude"]:.0f}°',
                showlegend=False,
                hovertemplate=(f'Start: {trace_data["longitude"]:.1f}°<br>'
                             f'R: {trace_data["radius"]:.3f} RE<br>'
                             f'Accuracy: ±{accuracy:.3f} RE<extra></extra>')
            ),
            row=1, col=1
        )
        
        # XY projection
        fig.add_trace(
            go.Scatter(
                x=x, y=y,
                mode='lines',
                line=dict(color=color, width=3),
                name=f'Closed {trace_data["longitude"]:.0f}°',
                showlegend=False,
                hovertemplate=f'Longitude: {trace_data["longitude"]:.1f}°<extra></extra>'
            ),
            row=1, col=2
        )
    
    # Plot open field lines
    for i, trace_data in enumerate(open_traces):
        if trace_data is None:
            continue
            
        trace = trace_data['trace']
        trace_coords = trace.get('combined', trace.get('forward', {})).get('coordinates', np.array([]))
        
        if len(trace_coords) < 2:
            continue
        
        x, y, z = trace_coords.T
        
        # Plot in red with dashed lines
        fig.add_trace(
            go.Scatter3d(
                x=x, y=y, z=z,
                mode='lines',
                line=dict(color='red', width=2, dash='dash'),
                name=f'Open {trace_data["longitude"]:.0f}°',
                showlegend=False,
                opacity=0.6
            ),
            row=1, col=1
        )
        
        fig.add_trace(
            go.Scatter(
                x=x, y=y,
                mode='lines',
                line=dict(color='red', width=2, dash='dash'),
                name=f'Open {trace_data["longitude"]:.0f}°',
                showlegend=False,
                opacity=0.6
            ),
            row=1, col=2
        )
    
    # Plot magnetopause boundary
    if show_magnetopause and len(magnetopause_points) > 0:
        mp_x, mp_y, mp_z = magnetopause_points.T
        
        # 3D scatter
        fig.add_trace(
            go.Scatter3d(
                x=mp_x, y=mp_y, z=mp_z,
                mode='markers',
                marker=dict(color='orange', size=10, symbol='square'),
                name='Magnetopause',
                showlegend=True
            ),
            row=1, col=1
        )
        
        # XY projection
        fig.add_trace(
            go.Scatter(
                x=mp_x, y=mp_y,
                mode='markers+lines',
                marker=dict(color='orange', size=10, symbol='square'),
                line=dict(color='orange', width=4),
                name='Magnetopause',
                showlegend=False
            ),
            row=1, col=2
        )
    
    # Plot radial distance vs longitude with error bars
    valid_mask = ~np.isnan(last_closed_r)
    if np.any(valid_mask):
        boundary_accuracy = results.get('boundary_accuracy', np.full_like(last_closed_r, 0.01))
        
        fig.add_trace(
            go.Scatter(
                x=longitudes[valid_mask], 
                y=last_closed_r[valid_mask],
                error_y=dict(
                    type='data',
                    array=boundary_accuracy[valid_mask],
                    visible=True,
                    thickness=2,
                    width=3
                ),
                mode='lines+markers',
                line=dict(color='blue', width=3),
                marker=dict(color='blue', size=8),
                name='Last Closed Radius',
                showlegend=True,
                hovertemplate=('Longitude: %{x:.1f}°<br>'
                             'Radius: %{y:.3f} RE<br>'
                             'Accuracy: ±%{error_y.array:.3f} RE<extra></extra>')
            ),
            row=2, col=2
        )
    
    # Add Earth
    # 3D Earth sphere
    u = np.linspace(0, 2 * np.pi, 20)
    v = np.linspace(0, np.pi, 20)
    x_earth = np.outer(np.cos(u), np.sin(v))
    y_earth = np.outer(np.sin(u), np.sin(v))
    z_earth = np.outer(np.ones(np.size(u)), np.cos(v))
    
    fig.add_trace(
        go.Surface(
            x=x_earth, y=y_earth, z=z_earth,
            colorscale=[[0, 'blue'], [1, 'blue']],
            opacity=0.3,
            showscale=False,
            name='Earth',
            hoverinfo='skip'
        ),
        row=1, col=1
    )
    
    # 2D Earth circle
    theta = np.linspace(0, 2*np.pi, 100)
    earth_x = np.cos(theta)
    earth_y = np.sin(theta)
    
    fig.add_trace(
        go.Scatter(
            x=earth_x, y=earth_y,
            mode='lines',
            line=dict(color='blue', width=4),
            name='Earth',
            hoverinfo='skip',
            showlegend=False
        ),
        row=1, col=2
    )
    
    # Update layout with aspect ratio 1
    coord_sys = results['summary']['coord_system']
    time_hours = results['summary']['time_hours']
    tolerance = results['summary'].get('tolerance', 0.01)
    
    fig.update_layout(
        title=f'Last Closed Field Lines (t={time_hours:.1f}h, tol={tolerance:.3f} RE)',
        height=800,
        width=1400,
        scene_aspectmode='data'  # This sets aspect ratio to 1 for 3D plot
    )
    
    # Update 3D scene
    fig.update_scenes(
        xaxis_title=f'X ({coord_sys}) [RE]',
        yaxis_title=f'Y ({coord_sys}) [RE]',
        zaxis_title=f'Z ({coord_sys}) [RE]'
    )
    
    # Update 2D plots with aspect ratio 1
    fig.update_xaxes(title_text=f'X ({coord_sys}) [RE]', row=1, col=2)
    fig.update_yaxes(title_text=f'Y ({coord_sys}) [RE]', row=1, col=2, 
                     scaleanchor="x", scaleratio=1)  # Aspect ratio 1
    
    fig.update_xaxes(title_text='Longitude [degrees]', row=2, col=2)
    fig.update_yaxes(title_text='Last Closed Radius [RE]', row=2, col=2)
    
    # Add grids
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
    
    return fig


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
    
    # Find last closed field lines
    results = find_last_closed_field_lines(
        ko, time_hours,
        n_traces=36,      # 36 traces at 10° intervals
        r_start=2.0,      # Start searching at 2 RE
        r_max=15.0,       # Search up to 15 RE
        r_step=0.5,       # 0.5 RE steps
        verbose=True
    )
    
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
    print(f"Found closed field lines: {summary['n_closed']}/{summary['n_total']}")
    if summary['n_closed'] > 0:
        print(f"Mean last closed radius: {summary['mean_radius']:.2f} RE")
        print(f"Range: [{summary['min_radius']:.2f}, {summary['max_radius']:.2f}] RE")
    
    return results


if __name__ == "__main__":
    print("=" * 60)
    print("LAST CLOSED MAGNETIC FIELD LINE FINDER")
    print("For use with Kamodo magnetospheric models")
    print("All coordinates in Earth radii (RE)")
    print("=" * 60)
    
    print("\n🚀 QUICK START GUIDE:")
    print("1. Load your Kamodo object: ko = your_kamodo_model")
    print("2. Find last closed field lines: results = find_last_closed_field_lines(ko, time_hours=12)")  
    print("3. Plot results: fig = plot_last_closed_field_lines(results)")
    
    print("\n🔧 USAGE EXAMPLE:")
    print("import kamodo_ccmc.tools.lastclosedmag as lcm")
    print("results = lcm.find_last_closed_field_lines(ko, time_hours=12.0)")
    print("lcm.plot_last_closed_field_lines(results).show()")

