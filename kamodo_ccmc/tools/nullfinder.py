"""
Magnetic Null Point Finder for Kamodo

Tools for finding and analyzing magnetic null points (reconnection sites) in magnetospheric models
All coordinates are in Earth radii (RE).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

# Import required tools from other modules
import kamodo_ccmc.tools.vectorfieldtracer as vft
import kamodo_ccmc.tools.lastclosedmag as lcm


def classify_field_line_topology(trace_result, earth_radius_threshold=1.05, 
                               imf_threshold=50.0, verbose=False):
    """
    Classify field line topology based on termination points and characteristics.
    
    Parameters:
    -----------
    trace_result : dict
        Result from trace_vector_line
    earth_radius_threshold : float
        Radius threshold for Earth surface detection (RE)
    imf_threshold : float
        Distance threshold for IMF classification (RE)
    verbose : bool
        Whether to print classification details
        
    Returns:
    --------
    topology : str
        'closed', 'open_north', 'open_south', 'imf_draped', 'unknown'
    details : dict
        Additional classification details
    """
    # Use combined trace if available
    if 'combined' in trace_result:
        coords = trace_result['combined']['coordinates']
        stop_reasons = trace_result['combined'].get('stop_reasons', [])
    else:
        # Combine forward and backward traces
        coords_list = []
        stop_reasons = []
        
        if 'backward' in trace_result:
            back_coords = trace_result['backward']['coordinates'][::-1]
            coords_list.extend(back_coords[:-1])  # Exclude duplicate middle point
            stop_reasons.append(trace_result['backward']['stop_reason'])
            
        if 'forward' in trace_result:
            forward_coords = trace_result['forward']['coordinates']
            coords_list.extend(forward_coords)
            stop_reasons.append(trace_result['forward']['stop_reason'])
            
        coords = np.array(coords_list) if coords_list else np.array([])
    
    if len(coords) < 2:
        return 'unknown', {'reason': 'insufficient_data'}
    
    # Get start and end points
    start_point = coords[0]
    end_point = coords[-1]
    
    # Calculate distances from Earth center
    start_r = np.linalg.norm(start_point)
    end_r = np.linalg.norm(end_point)
    
    # Check for Earth surface hits
    earth_hits = []
    if start_r < earth_radius_threshold:
        earth_hits.append('start')
    if end_r < earth_radius_threshold:
        earth_hits.append('end')
    
    # Count Earth surface hits from stop reasons
    earth_surface_hits = sum(1 for reason in stop_reasons if 'Earth_surface' in reason)
    
    details = {
        'start_point': start_point,
        'end_point': end_point,
        'start_r': start_r,
        'end_r': end_r,
        'stop_reasons': stop_reasons,
        'earth_surface_hits': earth_surface_hits,
        'max_distance': np.max(np.linalg.norm(coords, axis=1))
    }
    
    # Classification logic
    if earth_surface_hits >= 2:
        # Both ends hit Earth - closed field line
        topology = 'closed'
        details['classification_reason'] = 'both_ends_hit_earth'
        
    elif earth_surface_hits == 1:
        # One end hits Earth, other is open
        # Determine hemisphere based on the open end
        if len(coords) > 1:
            # Find which end doesn't hit Earth
            if start_r < earth_radius_threshold:
                open_end = end_point
            else:
                open_end = start_point
                
            # Classify based on Z component (GSM coordinates)
            if open_end[2] > 0:
                topology = 'open_north'
                details['classification_reason'] = 'open_to_north'
            else:
                topology = 'open_south'
                details['classification_reason'] = 'open_to_south'
        else:
            topology = 'unknown'
            details['classification_reason'] = 'insufficient_data'
            
    elif earth_surface_hits == 0:
        # Neither end hits Earth - could be IMF or draped
        max_r = np.max(np.linalg.norm(coords, axis=1))
        
        if max_r > imf_threshold:
            topology = 'imf_draped'
            details['classification_reason'] = 'far_from_earth'
        else:
            # Check field line curvature and direction
            # Simple heuristic: if field line is relatively straight, it's likely IMF
            if len(coords) > 5:
                # Calculate curvature
                path_length = np.sum(np.linalg.norm(np.diff(coords, axis=0), axis=1))
                direct_distance = np.linalg.norm(end_point - start_point)
                curvature_ratio = path_length / direct_distance if direct_distance > 0 else np.inf
                
                if curvature_ratio < 1.5:  # Relatively straight
                    topology = 'imf_draped'
                    details['classification_reason'] = 'straight_line_far_field'
                else:
                    topology = 'imf_draped'  # Curved but still far field
                    details['classification_reason'] = 'curved_far_field'
            else:
                topology = 'imf_draped'
                details['classification_reason'] = 'default_far_field'
    else:
        topology = 'unknown'
        details['classification_reason'] = 'unclear_termination'
    
    if verbose:
        print(f"Topology: {topology}")
        print(f"Reason: {details['classification_reason']}")
        print(f"Start: r={start_r:.2f} RE")
        print(f"End: r={end_r:.2f} RE")
        print(f"Stop reasons: {stop_reasons}")
    
    return topology, details


def find_precise_null_point_between_topologies(pos1, pos2, topo1, topo2, ko, time_hours, 
                                              tolerance=0.01, max_iterations=15, verbose=True):
    """
    Use bisection method to find precise null point location between two different topologies.
    
    Parameters:
    -----------
    pos1, pos2 : array_like
        Two positions with different topologies
    topo1, topo2 : str
        Topologies at pos1 and pos2
    ko : Kamodo object
        Kamodo object for field evaluation
    time_hours : float
        Time in hours
    tolerance : float
        Position tolerance in RE
    max_iterations : int
        Maximum bisection iterations
    verbose : bool
        Whether to print progress
        
    Returns:
    --------
    precise_null_point : dict
        Dictionary with precise null point information
    """
    if verbose:
        print(f"Finding precise null point between {topo1} and {topo2}")
        print(f"Position range: {np.linalg.norm(pos2 - pos1):.3f} RE")
    
    # Create tracer
    tracer = vft.KamodoVectorFieldTracer(
        ko, vector_components=('B_x', 'B_y', 'B_z'),
        field_type='magnetic'
    )
    
    # Bisection search
    p1, p2 = np.array(pos1), np.array(pos2)
    
    for iteration in range(max_iterations):
        # Midpoint
        p_mid = (p1 + p2) / 2
        distance = np.linalg.norm(p2 - p1)
        
        if distance < tolerance:
            if verbose:
                print(f"Converged after {iteration} iterations")
            break
        
        try:
            # Trace from midpoint
            trace_result = tracer.trace_vector_line(
                p_mid, time_hours,
                max_steps=3000,
                step_size_re=0.1,
                direction='both'
            )
            
            topology_mid, _ = classify_field_line_topology(trace_result, verbose=False)
            
            # Get field strength
            bx, by, bz = tracer.get_vector_field(p_mid[0], p_mid[1], p_mid[2], time_hours)
            if not np.isscalar(bx):
                bx, by, bz = float(bx[0]), float(by[0]), float(bz[0])
            field_strength = np.sqrt(bx**2 + by**2 + bz**2)
            
            if verbose and iteration < 5:
                print(f"  Iter {iteration}: topology={topology_mid}, "
                      f"field={field_strength:.2e} nT, dist={distance:.3f} RE")
            
            # Update search bounds
            if topology_mid == topo1:
                p1 = p_mid
            else:
                p2 = p_mid
                
        except Exception as e:
            if verbose:
                print(f"  Error at iteration {iteration}: {e}")
            break
    
    # Final evaluation at best position
    final_pos = (p1 + p2) / 2
    
    try:
        bx, by, bz = tracer.get_vector_field(final_pos[0], final_pos[1], final_pos[2], time_hours)
        if not np.isscalar(bx):
            bx, by, bz = float(bx[0]), float(by[0]), float(bz[0])
        final_field = np.sqrt(bx**2 + by**2 + bz**2)
        final_components = [bx, by, bz]
    except Exception:
        final_field = np.nan
        final_components = [np.nan, np.nan, np.nan]
    
    return {
        'location': final_pos,
        'field_strength': final_field,
        'field_components': final_components,
        'topologies': [topo1, topo2],
        'final_distance': np.linalg.norm(p2 - p1),
        'iterations': iteration + 1,
        'converged': distance < tolerance
    }


def find_magnetic_null_points(last_closed_results, ko, time_hours, 
                             radial_step=0.1, latitudinal_step=0.1,
                             max_lat_steps=20, max_radial_steps=5,
                             tolerance=0.05, verbose=True):
    """
    Find magnetic null points by searching along the magnetopause surface
    for topology transitions.
    
    Parameters:
    -----------
    last_closed_results : dict
        Results from find_last_closed_field_lines
    ko : Kamodo object
        Kamodo object for field evaluation
    time_hours : float
        Time in hours
    radial_step : float
        Initial radial step outward from last closed field line (RE)
    latitudinal_step : float
        Step size when moving North/South along magnetopause (RE)
    max_lat_steps : int
        Maximum steps North or South to search
    max_radial_steps : int
        Maximum radial steps outward to search
    tolerance : float
        Distance tolerance for bisection convergence (RE)
    verbose : bool
        Whether to print progress information
        
    Returns:
    --------
    null_point_candidates : list
        List of identified null point candidates
    topology_analysis : dict
        Detailed topology analysis results
    """
    if verbose:
        print("Finding magnetic null points along magnetopause surface...")
        print(f"Radial step: {radial_step} RE, Latitudinal step: {latitudinal_step} RE")
    
    # Get last closed field line data
    closed_traces = last_closed_results['closed_traces']
    
    if len(closed_traces) == 0:
        print("No closed field line traces found!")
        return [], {}
    
    # Create magnetic field tracer
    tracer = vft.KamodoVectorFieldTracer(
        ko, vector_components=('B_x', 'B_y', 'B_z'),
        field_type='magnetic'
    )
    
    # Find equatorial closed field lines (|Z| < 0.5 RE)
    equatorial_traces = []
    for trace_data in closed_traces:
        if trace_data is None:
            continue
        
        start_point = trace_data['start_point']
        if abs(start_point[2]) <= 0.5:  # Near equatorial plane
            equatorial_traces.append(trace_data)
    
    if len(equatorial_traces) == 0:
        print("No equatorial closed field lines found!")
        return [], {}
    
    if verbose:
        print(f"Found {len(equatorial_traces)} equatorial closed field lines")
    
    null_point_candidates = []
    detailed_analysis = []
    
    # Process each equatorial closed field line
    for i, trace_data in enumerate(equatorial_traces):
        longitude = trace_data['longitude']
        lon_rad = np.radians(longitude)
        base_radius = trace_data['radius']
        
        if verbose:
            print(f"\nProcessing longitude {longitude:.1f}° (trace {i+1}/{len(equatorial_traces)})")
            print(f"  Last closed radius: {base_radius:.2f} RE")
        
        # Step 1: Move outward radially from closed field line
        test_radius = base_radius + radial_step
        x_base = test_radius * np.cos(lon_rad)
        y_base = test_radius * np.sin(lon_rad)
        z_base = 0.0  # Start at equator
        
        try:
            # Test topology at equatorial point just outside magnetopause
            initial_trace = tracer.trace_vector_line(
                (x_base, y_base, z_base), time_hours,
                max_steps=3000, step_size_re=0.1, direction='both'
            )
            
            initial_topology, _ = classify_field_line_topology(initial_trace, verbose=False)
            
            if verbose:
                print(f"  Initial topology at ({x_base:.2f}, {y_base:.2f}, {z_base:.2f}): {initial_topology}")
            
            if initial_topology == 'closed':
                if verbose:
                    print(f"  Still closed at {test_radius:.2f} RE, skipping this longitude")
                continue
        except Exception as e:
            if verbose:
                print(f"  Error testing initial point: {e}")
            continue
        
        # Track if this is an IMF/draped field line initially
        found_imf_first = (initial_topology == 'imf_draped')
        
        if found_imf_first:
            if verbose:
                print(f"  Found IMF/draped field line initially. Searching for null point...")
                
            # If we found IMF first, we need to search North and South for open field lines
            search_directions = [1, -1]  # North and South
            open_points = {'open_north': None, 'open_south': None}
            
            for direction in search_directions:
                direction_name = "northward" if direction > 0 else "southward"
                
                if verbose:
                    print(f"    Searching {direction_name} for open field lines...")
                
                for step in range(1, max_lat_steps + 1):
                    # Move North or South from the IMF point
                    z_test = z_base + direction * step * latitudinal_step
                    
                    # Maintain constant radius from Earth center
                    r_meridional = np.sqrt(x_base**2 + z_test**2)
                    scale_factor = test_radius / r_meridional if r_meridional > 0 else 1.0
                    
                    x_test = x_base * scale_factor
                    y_test = y_base
                    z_test = z_test * scale_factor
                    
                    try:
                        test_trace = tracer.trace_vector_line(
                            (x_test, y_test, z_test), time_hours,
                            max_steps=3000, step_size_re=0.1, direction='both'
                        )
                        
                        test_topology, _ = classify_field_line_topology(test_trace, verbose=False)
                        
                        if test_topology in ['open_north', 'open_south']:
                            if verbose:
                                print(f"      Found {test_topology} at Z={z_test:.2f} RE")
                            
                            open_points[test_topology] = (x_test, y_test, z_test)
                            
                            # Stop once we find an open field line in this direction
                            break
                            
                    except Exception:
                        continue
            
            # Check if we found both open north and south
            if open_points['open_north'] is not None and open_points['open_south'] is not None:
                if verbose:
                    print("    Found both open_north and open_south - can locate null point")
                
                # Use bisection to find the null point between open_north and open_south
                north_pos = np.array(open_points['open_north'])
                south_pos = np.array(open_points['open_south'])
                
                precise_null = find_precise_null_point_between_topologies(
                    north_pos, south_pos, 'open_north', 'open_south',
                    ko, time_hours, tolerance=tolerance, verbose=verbose
                )
                
                null_point_candidates.append({
                    'location': precise_null['location'],
                    'field_strength': precise_null['field_strength'],
                    'field_components': precise_null['field_components'],
                    'topologies': ['open_north', 'open_south'],
                    'longitude': longitude,
                    'transition_type': 'open_north↔open_south',
                    'search_method': 'imf_to_open',
                    'quality_score': 1.0 / (precise_null['field_strength'] + 1e-12)
                })
                
                # Store the IMF position as well
                bx, by, bz = tracer.get_vector_field(x_base, y_base, z_base, time_hours)
                if not np.isscalar(bx):
                    bx, by, bz = float(bx[0]), float(by[0]), float(bz[0])
                imf_field_strength = np.sqrt(bx**2 + by**2 + bz**2)
                
                null_point_candidates.append({
                    'location': np.array([x_base, y_base, z_base]),
                    'field_strength': imf_field_strength,
                    'field_components': [bx, by, bz],
                    'topologies': ['imf_draped'],
                    'longitude': longitude,
                    'transition_type': 'imf_boundary',
                    'search_method': 'initial_imf',
                    'quality_score': 1.0 / (imf_field_strength + 1e-12)
                })
                
                continue  # Skip to next longitude
                
            elif verbose:
                print("    Could not find both open_north and open_south")
                
        # Normal case: open_north or open_south found first
        # Track search path
        topology_sequence = [{
            'position': np.array([x_base, y_base, z_base]),
            'topology': initial_topology,
            'step': 0
        }]
        
        # Determine search direction based on initial topology
        if initial_topology == 'open_north':
            # Search southward for open_south
            search_direction = -1
            target_topology = 'open_south'
            if verbose:
                print(f"  Found open_north, searching southward for open_south")
        elif initial_topology == 'open_south':
            # Search northward for open_north
            search_direction = 1
            target_topology = 'open_north'
            if verbose:
                print(f"  Found open_south, searching northward for open_north")
        else:
            # Already handled the IMF case, this must be 'unknown'
            if verbose:
                print(f"  Unhandled topology: {initial_topology}, skipping")
            continue
        
        # Search along the magnetopause surface
        found_transition = False
        transition_point = None
        
        for step in range(1, max_lat_steps + 1):
            # Move North or South along the magnetopause surface
            z_test = z_base + search_direction * step * latitudinal_step
            
            # Keep the same radial distance from Earth center
            r_meridional = np.sqrt(x_base**2 + z_test**2)
            scale_factor = test_radius / r_meridional if r_meridional > 0 else 1.0
            
            x_test = x_base * scale_factor
            y_test = y_base  # Y stays the same (moving in X-Z plane)
            z_test = z_test * scale_factor
            
            test_position = np.array([x_test, y_test, z_test])
            
            try:
                # Test topology at this position
                test_trace = tracer.trace_vector_line(
                    test_position, time_hours,
                    max_steps=3000, step_size_re=0.1, direction='both'
                )
                
                test_topology, _ = classify_field_line_topology(test_trace, verbose=False)
                
                # Get field strength
                bx, by, bz = tracer.get_vector_field(x_test, y_test, z_test, time_hours)
                if not np.isscalar(bx):
                    bx, by, bz = float(bx[0]), float(by[0]), float(bz[0])
                field_strength = np.sqrt(bx**2 + by**2 + bz**2)
                
                topology_sequence.append({
                    'position': test_position,
                    'topology': test_topology,
                    'step': step,
                    'field_strength': field_strength,
                    'field_components': [bx, by, bz]
                })
                
                if verbose and step <= 5:
                    print(f"    Step {step}: Z={z_test:.2f} RE, topology={test_topology}")
                
                # Check if we found the target topology
                if test_topology == target_topology:
                    found_transition = True
                    
                    # Found the transition! The null point is between the last two positions
                    pos1 = topology_sequence[-2]['position']  # Previous position
                    pos2 = topology_sequence[-1]['position']  # Current position
                    topo1 = topology_sequence[-2]['topology']
                    topo2 = topology_sequence[-1]['topology']
                    
                    # Use bisection to find precise transition point
                    precise_null = find_precise_null_point_between_topologies(
                        pos1, pos2, topo1, topo2, ko, time_hours,
                        tolerance=tolerance, verbose=verbose and i < 2
                    )
                    
                    if verbose:
                        print(f"  Found topology transition: {topo1} → {topo2}")
                        print(f"  Null point candidate: ({precise_null['location'][0]:.2f}, "
                              f"{precise_null['location'][1]:.2f}, {precise_null['location'][2]:.2f}) RE")
                        print(f"  Field strength: {precise_null['field_strength']:.2e} nT")
                    
                    null_point_candidates.append({
                        'location': precise_null['location'],
                        'field_strength': precise_null['field_strength'],
                        'field_components': precise_null['field_components'],
                        'topologies': [topo1, topo2],
                        'longitude': longitude,
                        'transition_type': f"{topo1}→{topo2}",
                        'search_direction': 'north' if search_direction > 0 else 'south',
                        'quality_score': 1.0 / (precise_null['field_strength'] + 1e-12),
                        'convergence_info': precise_null
                    })
                    
                    transition_point = precise_null['location']
                    break
                    
            except Exception as e:
                if verbose:
                    print(f"    Error at step {step}: {e}")
                continue
        
        if not found_transition:
            if verbose:
                print(f"  No topology transition found within {max_lat_steps} steps")
            continue
        
        # Search for IMF/draped field lines farther out
        if transition_point is not None:
            if verbose:
                print(f"  Searching for IMF/draped field lines farther out...")
            
            # Use the null point as reference
            null_radius = np.linalg.norm(transition_point)
            null_lon = np.arctan2(transition_point[1], transition_point[0])
            null_lat = np.arcsin(transition_point[2] / null_radius) if null_radius > 0 else 0
            
            for radial_step_out in range(1, max_radial_steps + 1):
                imf_radius = null_radius + radial_step_out * radial_step
                
                # Keep the same angular position as the null point
                x_imf = imf_radius * np.cos(null_lat) * np.cos(null_lon)
                y_imf = imf_radius * np.cos(null_lat) * np.sin(null_lon)  
                z_imf = imf_radius * np.sin(null_lat)
                
                try:
                    imf_trace = tracer.trace_vector_line(
                        (x_imf, y_imf, z_imf), time_hours,
                        max_steps=2000, step_size_re=0.15, direction='both'
                    )
                    
                    imf_topology, _ = classify_field_line_topology(imf_trace, verbose=False)
                    
                    if imf_topology == 'imf_draped':
                        # Get field strength
                        bx, by, bz = tracer.get_vector_field(x_imf, y_imf, z_imf, time_hours)
                        if not np.isscalar(bx):
                            bx, by, bz = float(bx[0]), float(by[0]), float(bz[0])
                        imf_field_strength = np.sqrt(bx**2 + by**2 + bz**2)
                        
                        null_point_candidates.append({
                            'location': np.array([x_imf, y_imf, z_imf]),
                            'field_strength': imf_field_strength,
                            'field_components': [bx, by, bz],
                            'topologies': ['imf_draped'],
                            'longitude': longitude,
                            'transition_type': 'imf_boundary',
                            'quality_score': 1.0 / (imf_field_strength + 1e-12),
                            'radial_distance_from_null': radial_step_out * radial_step
                        })
                        
                        if verbose:
                            print(f"    Found IMF/draped at R={imf_radius:.2f} RE, "
                                  f"field={imf_field_strength:.2e} nT")
                        break
                        
                except Exception as e:
                    if verbose and radial_step_out == 1:
                        print(f"    Error searching for IMF: {e}")
                    continue
        
        # Store detailed analysis for this longitude
        detailed_analysis.append({
            'longitude': longitude,
            'base_radius': base_radius,
            'initial_topology': initial_topology,
            'topology_sequence': topology_sequence,
            'found_transition': found_transition,
            'transition_point': transition_point
        })
    
    # Sort candidates by quality (lowest field strength first)
    null_point_candidates.sort(key=lambda x: x['field_strength'])
    
    # Create topology analysis summary
    topology_analysis = {
        'detailed_analysis': detailed_analysis,
        'n_candidates': len(null_point_candidates),
        'search_parameters': {
            'radial_step': radial_step,
            'latitudinal_step': latitudinal_step,
            'max_lat_steps': max_lat_steps,
            'max_radial_steps': max_radial_steps,
            'tolerance': tolerance
        }
    }
    
    if verbose:
        print(f"\nNull point search results:")
        print(f"  Processed {len(equatorial_traces)} equatorial field lines")
        print(f"  Found topology transitions: {sum(1 for a in detailed_analysis if a['found_transition'])}")
        print(f"  Total null point candidates: {len(null_point_candidates)}")
        
        # Categorize candidates
        topology_transitions = [c for c in null_point_candidates if '→' in c['transition_type']]
        imf_boundaries = [c for c in null_point_candidates if c['transition_type'] == 'imf_boundary']
        
        print(f"  Topology transition null points: {len(topology_transitions)}")
        print(f"  IMF boundary points: {len(imf_boundaries)}")
        
        if null_point_candidates:
            best = null_point_candidates[0]
            print(f"\nBest null point candidate:")
            print(f"  Location: ({best['location'][0]:.2f}, "
                  f"{best['location'][1]:.2f}, {best['location'][2]:.2f}) RE")
            print(f"  Field strength: {best['field_strength']:.2e} nT")
            print(f"  Transition: {best['transition_type']}")
            print(f"  Longitude: {best['longitude']:.1f}°")
    
    return null_point_candidates, topology_analysis


def refine_null_point_location(null_point_candidate, ko, time_hours, search_radius=0.5, 
                              n_steps=10, verbose=True):
    """
    Refine the null point location by finding the minimum field strength in the vicinity.
    
    Parameters:
    -----------
    null_point_candidate : dict
        Initial null point candidate from find_magnetic_null_points
    ko : Kamodo object
        Kamodo object for field evaluation
    time_hours : float
        Time in hours
    search_radius : float
        Search radius around initial point (RE)
    n_steps : int
        Number of refinement steps
    verbose : bool
        Whether to print progress information
        
    Returns:
    --------
    refined_null_point : dict
        Refined null point with updated location and field values
    """
    initial_location = null_point_candidate['location']
    
    if verbose:
        print(f"Refining null point location starting from: "
              f"({initial_location[0]:.2f}, {initial_location[1]:.2f}, {initial_location[2]:.2f}) RE")
        print(f"Initial field strength: {null_point_candidate['field_strength']:.2e} nT")
    
    # Create tracer for field evaluation
    tracer = vft.KamodoVectorFieldTracer(
        ko, vector_components=('B_x', 'B_y', 'B_z'),
        field_type='magnetic'
    )
    
    # Function to evaluate field strength at a point
    def field_strength(position):
        try:
            bx, by, bz = tracer.get_vector_field(position[0], position[1], position[2], time_hours)
            
            # Handle array returns
            if not np.isscalar(bx):
                bx, by, bz = float(bx[0]), float(by[0]), float(bz[0])
            
            return np.sqrt(bx**2 + by**2 + bz**2)
        except Exception:
            return 1e10  # Return very large value if evaluation fails
    
    # Objective function for minimization
    def objective_function(position):
        return field_strength(position)
    
    # Set bounds for optimization
    bounds = [(initial_location[i] - search_radius, initial_location[i] + search_radius) 
              for i in range(3)]
    
    # Multi-start optimization for robustness
    best_position = initial_location
    best_field = field_strength(initial_location)
    
    # First do a grid search to find good starting points
    n_grid = 5
    grid_points = []
    for i in range(n_grid):
        x = initial_location[0] - search_radius + 2 * search_radius * i / (n_grid - 1)
        for j in range(n_grid):
            y = initial_location[1] - search_radius + 2 * search_radius * j / (n_grid - 1)
            for k in range(n_grid):
                z = initial_location[2] - search_radius + 2 * search_radius * k / (n_grid - 1)
                grid_points.append(np.array([x, y, z]))
    
    # Evaluate field at grid points
    grid_results = []
    for point in grid_points:
        field = field_strength(point)
        grid_results.append((point, field))
    
    # Sort by field strength
    grid_results.sort(key=lambda x: x[1])
    
    # Use best grid points as starting points for optimization
    start_points = [res[0] for res in grid_results[:5]]
    
    # Run optimization from each starting point
    for i, start_pos in enumerate(start_points):
        try:
            result = minimize(
                objective_function, 
                start_pos, 
                method='L-BFGS-B',
                bounds=bounds,
                options={'maxiter': 50, 'ftol': 1e-10}
            )
            
            if result.success:
                final_position = result.x
                final_field = result.fun
                
                if final_field < best_field:
                    best_position = final_position
                    best_field = final_field
                    
                    if verbose:
                        print(f"  Improved solution: ({final_position[0]:.3f}, "
                              f"{final_position[1]:.3f}, {final_position[2]:.3f}) RE, "
                              f"field = {final_field:.2e} nT")
            
        except Exception as e:
            if verbose:
                print(f"  Optimization attempt {i+1} failed: {e}")
    
    # Evaluate field components at best position
    try:
        bx, by, bz = tracer.get_vector_field(
            best_position[0], best_position[1], best_position[2], time_hours)
        
        # Handle array returns
        if not np.isscalar(bx):
            bx, by, bz = float(bx[0]), float(by[0]), float(bz[0])
            
    except Exception:
        bx, by, bz = np.nan, np.nan, np.nan
    
    # Create refined null point
    refined_null_point = {
        'location': best_position,
        'field_strength': best_field,
        'field_components': [bx, by, bz],
        'topologies': null_point_candidate['topologies'],
        'original_location': initial_location,
        'improvement': null_point_candidate['field_strength'] / best_field,
        'search_radius': search_radius
    }
    
    if verbose:
        print(f"Refinement complete:")
        print(f"  Original: ({initial_location[0]:.3f}, "
              f"{initial_location[1]:.3f}, {initial_location[2]:.3f}) RE, "
              f"field = {null_point_candidate['field_strength']:.2e} nT")
        print(f"  Refined: ({best_position[0]:.3f}, "
              f"{best_position[1]:.3f}, {best_position[2]:.3f}) RE, "
              f"field = {best_field:.2e} nT")
        print(f"  Improvement factor: {refined_null_point['improvement']:.2f}x")
    
    return refined_null_point


def trace_special_field_lines_near_null(null_point, ko, time_hours, 
                                       n_lines=12, radius=0.5, max_steps=5000):
    """
    Trace field lines in the vicinity of the null point to visualize its structure.
    
    Parameters:
    -----------
    null_point : dict
        Null point information
    ko : Kamodo object
        Kamodo object for field evaluation
    time_hours : float
        Time in hours
    n_lines : int
        Number of field lines to trace
    radius : float
        Radius around null point for starting points (RE)
    max_steps : int
        Maximum integration steps
    
    Returns:
    --------
    special_traces : dict
        Dictionary of traces around the null point
    """
    null_location = null_point['location']
    
    # Create magnetic field tracer
    tracer = vft.KamodoVectorFieldTracer(
        ko, vector_components=('B_x', 'B_y', 'B_z'),
        field_type='magnetic'
    )
    
    # Create starting points around null point in a spherical distribution
    start_points = []
    
    # Use golden spiral distribution for even coverage
    golden_angle = np.pi * (3 - np.sqrt(5))  # Golden angle in radians
    
    for i in range(n_lines):
        # Golden spiral distribution
        y = 1 - (i / (n_lines - 1)) * 2 if n_lines > 1 else 0  # y goes from 1 to -1
        radius_at_y = np.sqrt(1 - y**2)  # radius at y
        theta = golden_angle * i  # golden angle increment
        
        x = np.cos(theta) * radius_at_y
        z = np.sin(theta) * radius_at_y
        
        # Scale and offset
        point = null_location + radius * np.array([x, y, z])
        start_points.append(point)
    
    # Trace field lines from these starting points
    special_traces = []
    
    for i, start_point in enumerate(start_points):
        try:
            trace_result = tracer.trace_vector_line(
                start_point, time_hours,
                max_steps=max_steps,
                step_size_re=0.05,  # Smaller step size for detail
                direction='both',
                adaptive_step=True
            )
            
            # Classify topology
            topology, details = classify_field_line_topology(trace_result)
            
            special_traces.append({
                'trace': trace_result,
                'start_point': start_point,
                'distance_from_null': np.linalg.norm(start_point - null_location),
                'topology': topology,
                'details': details
            })
            
        except Exception as e:
            print(f"Error tracing from point {i+1}: {e}")
    
    # Group traces by topology
    topology_groups = {}
    for trace_info in special_traces:
        topology = trace_info['topology']
        if topology not in topology_groups:
            topology_groups[topology] = []
        topology_groups[topology].append(trace_info)
    
    return {
        'special_traces': special_traces,
        'topology_groups': topology_groups,
        'null_point': null_point,
        'n_lines': n_lines,
        'radius': radius
    }


def plot_null_point_structure(null_structure, backend='plotly', show_separatrix=True):
    """
    Plot the magnetic field structure around a null point.
    
    Parameters:
    -----------
    null_structure : dict
        Output from trace_special_field_lines_near_null
    backend : str
        'matplotlib' or 'plotly'
    show_separatrix : bool
        Whether to highlight the separatrix surfaces
        
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
            return _plot_null_structure_plotly(null_structure, show_separatrix)
    
    return _plot_null_structure_matplotlib(null_structure, show_separatrix)


def _plot_null_structure_matplotlib(null_structure, show_separatrix=True):
    """Create matplotlib plot of null point magnetic field structure."""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    special_traces = null_structure['special_traces']
    topology_groups = null_structure['topology_groups']
    null_point = null_structure['null_point']
    
    # Color map for topologies
    topology_colors = {
        'closed': 'blue',
        'open_north': 'red',
        'open_south': 'green',
        'imf_draped': 'orange',
        'unknown': 'gray',
        'error': 'black'
    }
    
    # Plot field lines by topology group
    for topology, traces in topology_groups.items():
        color = topology_colors.get(topology, 'gray')
        
        for i, trace_info in enumerate(traces):
            trace_result = trace_info['trace']
            
            # Use combined trace if available
            if 'combined' in trace_result:
                coords = trace_result['combined']['coordinates']
            elif 'forward' in trace_result:
                coords = trace_result['forward']['coordinates']
            else:
                coords = trace_result['backward']['coordinates']
            
            if len(coords) < 2:
                continue
            
            x, y, z = coords.T
            
            ax.plot(x, y, z, c=color, linewidth=1.5, alpha=0.7,
                   label=topology if i == 0 else None)
    
    # Plot null point
    null_location = null_point['location']
    ax.scatter([null_location[0]], [null_location[1]], [null_location[2]],
               c='purple', s=100, marker='D', alpha=1.0,
               edgecolors='white', linewidth=1, label='Null Point')
    
    # Add Earth
    u = np.linspace(0, 2 * np.pi, 20)
    v = np.linspace(0, np.pi, 20)
    x_earth = np.outer(np.cos(u), np.sin(v))
    y_earth = np.outer(np.sin(u), np.sin(v))
    z_earth = np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x_earth, y_earth, z_earth, color='lightblue', alpha=0.3)
    
    # Set equal aspect ratio
    all_coords = []
    for trace_info in special_traces:
        trace_result = trace_info['trace']
        if 'combined' in trace_result:
            all_coords.extend(trace_result['combined']['coordinates'])
        elif 'forward' in trace_result:
            all_coords.extend(trace_result['forward']['coordinates'])
        elif 'backward' in trace_result:
            all_coords.extend(trace_result['backward']['coordinates'])
    
    if all_coords:
        all_coords = np.array(all_coords)
        max_range = np.array([all_coords[:, 0].max()-all_coords[:, 0].min(),
                             all_coords[:, 1].max()-all_coords[:, 1].min(),
                             all_coords[:, 2].max()-all_coords[:, 2].min()]).max() / 2.0
        mid_x = (all_coords[:, 0].max()+all_coords[:, 0].min()) * 0.5
        mid_y = (all_coords[:, 1].max()+all_coords[:, 1].min()) * 0.5
        mid_z = (all_coords[:, 2].max()+all_coords[:, 2].min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    # Set labels and title
    ax.set_xlabel('X [RE]')
    ax.set_ylabel('Y [RE]')
    ax.set_zlabel('Z [RE]')
    ax.set_title('Magnetic Field Structure Around Null Point')
    
    # Add legend with unique entries
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='upper right')
    
    plt.tight_layout()
    return fig


def _plot_null_structure_plotly(null_structure, show_separatrix=True):
    """Create Plotly plot of null point magnetic field structure."""
    
    special_traces = null_structure['special_traces']
    topology_groups = null_structure['topology_groups']
    null_point = null_structure['null_point']
    
    # Create plot
    fig = make_subplots(
        rows=1, cols=1,
        specs=[[{"type": "scene"}]],
        subplot_titles=["Magnetic Field Structure Around Null Point"]
    )
    
    # Color map for topologies
    topology_colors = {
        'closed': 'blue',
        'open_north': 'red',
        'open_south': 'green',
        'imf_draped': 'orange',
        'unknown': 'gray',
        'error': 'black'
    }
    
    # Plot field lines by topology group
    for topology, traces in topology_groups.items():
        color = topology_colors.get(topology, 'gray')
        
        for i, trace_info in enumerate(traces):
            trace_result = trace_info['trace']
            
            # Use combined trace if available
            if 'combined' in trace_result:
                coords = trace_result['combined']['coordinates']
            elif 'forward' in trace_result:
                coords = trace_result['forward']['coordinates']
            else:
                coords = trace_result['backward']['coordinates']
            
            if len(coords) < 2:
                continue
            
            x, y, z = coords.T
            
            fig.add_trace(
                go.Scatter3d(
                    x=x, y=y, z=z,
                    mode='lines',
                    line=dict(
                        color=color,
                        width=3,
                        dash='solid'
                    ),
                    opacity=0.7,
                    name=f"{topology.replace('_', ' ').title()}",
                    showlegend=(i == 0),
                    hovertemplate=f"{topology.replace('_', ' ').title()}<extra></extra>"
                )
            )
    
    # Plot null point
    null_location = null_point['location']
    fig.add_trace(
        go.Scatter3d(
            x=[null_location[0]],
            y=[null_location[1]],
            z=[null_location[2]],
            mode='markers',
            marker=dict(
                color='purple',
                size=10,
                symbol='diamond',
                line=dict(color='white', width=2)
            ),
            name="Null Point",
            showlegend=True,
            hovertemplate=(f"Null Point<br>"
                         f"X: {null_location[0]:.3f} RE<br>"
                         f"Y: {null_location[1]:.3f} RE<br>"
                         f"Z: {null_location[2]:.3f} RE<br>"
                         f"Field: {null_point['field_strength']:.2e} nT<extra></extra>")
        )
    )
    
    # Plot Earth
    u = np.linspace(0, 2 * np.pi, 20)
    v = np.linspace(0, np.pi, 20)
    x_earth = np.outer(np.cos(u), np.sin(v))
    y_earth = np.outer(np.sin(u), np.sin(v))
    z_earth = np.outer(np.ones(np.size(u)), np.cos(v))
    
    fig.add_trace(
        go.Surface(
            x=x_earth, y=y_earth, z=z_earth,
            colorscale=[[0, 'lightblue'], [1, 'lightblue']],
            opacity=0.3,
            showscale=False,
            name='Earth',
            hoverinfo='skip'
        )
    )
    
    # Update layout
    fig.update_layout(
        title=f"Magnetic Field Structure Around Null Point",
        scene=dict(
            xaxis_title="X [RE]",
            yaxis_title="Y [RE]",
            zaxis_title="Z [RE]",
            aspectmode='data'
        ),
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        ),
        margin=dict(r=20, l=10, b=10, t=50),
        height=800,
        width=1000
    )
    
    return fig


def plot_null_point_analysis(null_point_candidates, topology_analysis, last_closed_results, 
                            backend='plotly', n_best=5):
    """
    Plot the null point analysis results.
    
    Parameters:
    -----------
    null_point_candidates : list
        List of null point candidates
    topology_analysis : dict
        Topology analysis results
    last_closed_results : dict
        Results from last closed field line analysis
    backend : str
        'matplotlib' or 'plotly'
    n_best : int
        Number of best null points to highlight
        
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
            return _plot_null_analysis_plotly(null_point_candidates, topology_analysis, 
                                            last_closed_results, n_best)
    
    return _plot_null_analysis_matplotlib(null_point_candidates, topology_analysis, 
                                        last_closed_results, n_best)


def _plot_null_analysis_matplotlib(null_point_candidates, topology_analysis, 
                                 last_closed_results, n_best=5):
    """Create matplotlib plot of null point analysis."""
    fig = plt.figure(figsize=(16, 12))
    
    # Create subplots
    ax1 = fig.add_subplot(221, projection='3d')  # 3D view
    ax2 = fig.add_subplot(222)  # XY topology map
    ax3 = fig.add_subplot(223)  # Field strength
    ax4 = fig.add_subplot(224)  # Statistics
    
    # Plot last closed field lines
    closed_traces = last_closed_results['closed_traces']
    for trace_data in closed_traces[:10]:  # Limit for clarity
        if trace_data is None:
            continue
            
        trace = trace_data['trace']
        coords = trace.get('combined', trace.get('forward', {})).get('coordinates', np.array([]))
        
        if len(coords) < 2:
            continue
        
        x, y, z = coords.T
        ax1.plot(x, y, z, color='lightblue', alpha=0.5, linewidth=1.0)
        ax2.plot(x, y, color='lightblue', alpha=0.5, linewidth=1.0)
    
    # Plot null point candidates
    if null_point_candidates:
        null_locations = np.array([np['location'] for np in null_point_candidates[:n_best]])
        null_strengths = [np['field_strength'] for np in null_point_candidates[:n_best]]
        
        # Size markers by inverse field strength
        max_strength = max(null_strengths) if null_strengths else 1
        marker_sizes = [200 * (max_strength / (strength + 1e-12)) for strength in null_strengths]
        marker_sizes = np.clip(marker_sizes, 50, 400)
        
        # 3D null points
        ax1.scatter(null_locations[:, 0], null_locations[:, 1], null_locations[:, 2],
                   c='purple', s=marker_sizes, marker='D', alpha=1.0,
                   edgecolors='white', linewidth=1, label='Null Points')
        
        # 2D null points
        ax2.scatter(null_locations[:, 0], null_locations[:, 1],
                   c='purple', s=marker_sizes, marker='D', alpha=1.0,
                   edgecolors='white', linewidth=1)
        
        # Field strength plot
        ax3.semilogy(range(1, len(null_strengths) + 1), null_strengths, 'o-', c='purple')
        ax3.set_xlabel('Null Point Candidate Rank')
        ax3.set_ylabel('Magnetic Field Strength [nT]')
        ax3.set_title('Field Strength at Null Point Candidates')
        ax3.grid(True, alpha=0.3)
    
    # Statistics
    if null_point_candidates:
        transition_types = [c.get('transition_type', 'unknown') for c in null_point_candidates]
        unique_types, counts = np.unique(transition_types, return_counts=True)
        
        ax4.pie(counts, labels=unique_types, autopct='%1.1f%%', startangle=90)
        ax4.set_title('Null Point Types')
    
    # Add Earth
    u = np.linspace(0, 2 * np.pi, 20)
    v = np.linspace(0, np.pi, 20)
    x_earth = np.outer(np.cos(u), np.sin(v))
    y_earth = np.outer(np.sin(u), np.sin(v))
    z_earth = np.outer(np.ones(np.size(u)), np.cos(v))
    ax1.plot_surface(x_earth, y_earth, z_earth, color='lightblue', alpha=0.3)
    
    # 2D Earth circle
    theta = np.linspace(0, 2*np.pi, 100)
    ax2.plot(np.cos(theta), np.sin(theta), c='lightblue', linewidth=2)
    ax2.set_aspect('equal')
    
    # Set labels
    coord_sys = last_closed_results['summary']['coord_system']
    time_hours = last_closed_results['summary']['time_hours']
    
    ax1.set_xlabel(f'X ({coord_sys}) [RE]')
    ax1.set_ylabel(f'Y ({coord_sys}) [RE]')
    ax1.set_zlabel(f'Z ({coord_sys}) [RE]')
    ax1.set_title(f'Null Point Analysis 3D (t={time_hours:.1f}h)')
    
    ax2.set_xlabel(f'X ({coord_sys}) [RE]')
    ax2.set_ylabel(f'Y ({coord_sys}) [RE]')
    ax2.set_title('Equatorial Plane View')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig


def _plot_null_analysis_plotly(null_point_candidates, topology_analysis, 
                             last_closed_results, n_best=5):
    """Create Plotly plot of null point analysis."""
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{"type": "scene", "rowspan": 2}, {"type": "xy"}],
               [None, {"type": "xy"}]],
        subplot_titles=[
            'Null Point Analysis 3D',
            'Equatorial Plane View',
            'Field Strength at Null Points'
        ],
        horizontal_spacing=0.1,
        vertical_spacing=0.1
    )
    
    # Plot last closed field lines
    closed_traces = last_closed_results['closed_traces']
    for i, trace_data in enumerate(closed_traces[:10]):  # Limit for clarity
        if trace_data is None:
            continue
            
        trace = trace_data['trace']
        coords = trace.get('combined', trace.get('forward', {})).get('coordinates', np.array([]))
        
        if len(coords) < 2:
            continue
        
        x, y, z = coords.T
        
        fig.add_trace(
            go.Scatter3d(
                x=x, y=y, z=z,
                mode='lines',
                line=dict(color='lightblue', width=2),
                name='Last Closed Field Lines',
                showlegend=i == 0,
                opacity=0.6
            ),
            row=1, col=1
        )
        
        fig.add_trace(
            go.Scatter(
                x=x, y=y,
                mode='lines',
                line=dict(color='lightblue', width=2),
                name='Last Closed (2D)',
                showlegend=False,
                opacity=0.6
            ),
            row=1, col=2
        )
    
    # Plot null point candidates
    if null_point_candidates:
        null_locations = np.array([np['location'] for np in null_point_candidates[:n_best]])
        null_strengths = [np['field_strength'] for np in null_point_candidates[:n_best]]
        
        # Size markers by inverse field strength
        max_strength = max(null_strengths) if null_strengths else 1
        marker_sizes = [20 * (max_strength / (strength + 1e-12)) for strength in null_strengths]
        marker_sizes = np.clip(marker_sizes, 5, 30)
        
        # 3D null points
        fig.add_trace(
            go.Scatter3d(
                x=null_locations[:, 0], y=null_locations[:, 1], z=null_locations[:, 2],
                mode='markers',
                marker=dict(
                    color='purple',
                    size=marker_sizes,
                    symbol='diamond',
                    line=dict(color='white', width=2)
                ),
                name='Null Point Candidates',
                showlegend=True,
                hovertemplate=('Null Point<br>'
                             'X: %{x:.2f} RE<br>'
                             'Y: %{y:.2f} RE<br>'
                             'Z: %{z:.2f} RE<br>'
                             'Field: %{text}<extra></extra>'),
                text=[f'{strength:.2e} nT' for strength in null_strengths]
            ),
            row=1, col=1
        )
        
        # 2D null points
        fig.add_trace(
            go.Scatter(
                x=null_locations[:, 0], y=null_locations[:, 1],
                mode='markers',
                marker=dict(
                    color='purple',
                    size=marker_sizes,
                    symbol='diamond',
                    line=dict(color='white', width=2)
                ),
                name='Null Points (2D)',
                showlegend=False,
                hovertemplate=('Null Point<br>'
                             'X: %{x:.2f} RE<br>'
                             'Y: %{y:.2f} RE<br>'
                             'Field: %{text}<extra></extra>'),
                text=[f'{strength:.2e} nT' for strength in null_strengths]
            ),
            row=1, col=2
        )
        
        # Field strength plot
        fig.add_trace(
            go.Scatter(
                x=list(range(1, len(null_strengths) + 1)),
                y=null_strengths,
                mode='markers+lines',
                marker=dict(color='purple', size=8),
                line=dict(color='purple', width=2),
                name='Field Strength',
                showlegend=False
            ),
            row=2, col=2
        )
    
    # Add Earth
    u = np.linspace(0, 2 * np.pi, 20)
    v = np.linspace(0, np.pi, 20)
    x_earth = np.outer(np.cos(u), np.sin(v))
    y_earth = np.outer(np.sin(u), np.sin(v))
    z_earth = np.outer(np.ones(np.size(u)), np.cos(v))
    
    fig.add_trace(
        go.Surface(
            x=x_earth, y=y_earth, z=z_earth,
            colorscale=[[0, 'lightblue'], [1, 'lightblue']],
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
            line=dict(color='lightblue', width=4),
            name='Earth',
            hoverinfo='skip',
            showlegend=False
        ),
        row=1, col=2
    )
    
    # Update layout
    coord_sys = last_closed_results['summary']['coord_system']
    time_hours = last_closed_results['summary']['time_hours']
    
    fig.update_layout(
        title=f'Null Point Analysis (t={time_hours:.1f}h)',
        height=800,
        width=1400,
        scene_aspectmode='data'
    )
    
    # Update 3D scene
    fig.update_scenes(
        xaxis_title=f'X ({coord_sys}) [RE]',
        yaxis_title=f'Y ({coord_sys}) [RE]',
        zaxis_title=f'Z ({coord_sys}) [RE]'
    )
    
    # Update 2D plots
    fig.update_xaxes(title_text=f'X ({coord_sys}) [RE]', row=1, col=2)
    fig.update_yaxes(title_text=f'Y ({coord_sys}) [RE]', row=1, col=2, 
                     scaleanchor="x", scaleratio=1)
    
    fig.update_xaxes(title_text='Null Point Rank', row=2, col=2)
    fig.update_yaxes(title_text='Field Strength [nT]', row=2, col=2, type='log')
    
    # Add grids
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
    
    return fig


def example_find_magnetic_nulls(ko, time_hours=12.0):
    """
    Complete example of finding magnetic null points.
    
    Parameters:
    -----------
    ko : Kamodo object
        Your Kamodo object
    time_hours : float
        Time in hours
    
    Returns:
    --------
    results : dict
        Dictionary containing analysis results
    """
    print("Example: Finding Magnetic Null Points")
    print("=" * 50)
    
    # Step 1: Find last closed field lines
    print("\nStep 1: Finding last closed field lines...")
    
    last_closed_results = lcm.find_last_closed_field_lines(
        ko, time_hours,
        n_traces=36,
        r_start=2.0,
        r_max=15.0,
        tolerance=0.02,
        verbose=True
    )
    
    # Step 2: Find magnetic null points
    print("\nStep 2: Finding magnetic null points...")
    null_point_candidates, topology_analysis = find_magnetic_null_points(
        last_closed_results, ko, time_hours,
        radial_step=0.1,
        latitudinal_step=0.1,
        max_lat_steps=20,
        max_radial_steps=5,
        verbose=True
    )
    
    if not null_point_candidates:
        print("No null point candidates found.")
        return {
            'last_closed_results': last_closed_results,
            'null_point_candidates': [],
            'topology_analysis': topology_analysis
        }
    
    # Step 3: Refine the best candidates
    print("\nStep 3: Refining null point locations...")
    refined_candidates = []
    
    for i, candidate in enumerate(null_point_candidates[:3]):  # Top 3 candidates
        print(f"\nRefining candidate {i+1}:")
        
        refined = refine_null_point_location(
            candidate, ko, time_hours,
            search_radius=0.2,
            verbose=True
        )
        refined_candidates.append(refined)
    
    # Step 4: Analyze field structure around best null point
    print("\nStep 4: Analyzing field structure around best null point...")
    best_null = refined_candidates[0]
    
    null_structure = trace_special_field_lines_near_null(
        best_null, ko, time_hours,
        n_lines=24,
        radius=0.3
    )
    
    # Step 5: Create visualizations
    print("\nStep 5: Creating visualizations...")
    
    # Analysis plot
    analysis_fig = plot_null_point_analysis(
        null_point_candidates, topology_analysis, last_closed_results,
        backend='plotly'
    )
    
    # Structure plot
    structure_fig = plot_null_point_structure(
        null_structure, backend='plotly'
    )
    
    # Show plots if in notebook environment
    try:
        from IPython import get_ipython
        if get_ipython() is not None:
            print("\nDisplaying interactive plots...")
            analysis_fig.show()
            structure_fig.show()
    except ImportError:
        print("\nTo view interactive plots, call:")
        print("  analysis_fig.show()")
        print("  structure_fig.show()")
    
    # Summary
    print(f"\nSummary:")
    print(f"  Magnetic null points found: {len(null_point_candidates)}")
    print(f"  Best null point location: ({best_null['location'][0]:.2f}, "
          f"{best_null['location'][1]:.2f}, {best_null['location'][2]:.2f}) RE")
    print(f"  Field strength at null: {best_null['field_strength']:.2e} nT")
    print(f"  Improvement from refinement: {best_null['improvement']:.1f}x")
    
    # Return all results
    return {
        'last_closed_results': last_closed_results,
        'null_point_candidates': null_point_candidates,
        'topology_analysis': topology_analysis,
        'refined_candidates': refined_candidates,
        'best_null_point': best_null,
        'null_structure': null_structure,
        'analysis_fig': analysis_fig,
        'structure_fig': structure_fig
    }


if __name__ == "__main__":
    print("=" * 60)
    print("MAGNETIC NULL POINT FINDER")
    print("For use with Kamodo magnetospheric models")
    print("All coordinates in Earth radii (RE)")
    print("=" * 60)
    
    print("\n🚀 QUICK START GUIDE:")
    print("1. Load your Kamodo object: ko = your_kamodo_model")
    print("2. Find last closed field lines: results = lcm.find_last_closed_field_lines(ko, time_hours)")
    print("3. Find null points: nulls, analysis = find_magnetic_null_points(results, ko, time_hours)")
    print("4. Plot results: fig = plot_null_point_analysis(nulls, analysis, results)")
    
    print("\n🔧 USAGE EXAMPLE:")
    print("import kamodo_ccmc.tools.nullfinder as nf")
    print("import kamodo_ccmc.tools.lastclosedmag as lcm")
    print("results = nf.example_find_magnetic_nulls(ko, time_hours=12.0)")
    print("results['analysis_fig'].show()")
    print("results['structure_fig'].show()")
    
    print("\n📊 KEY FEATURES:")
    print("• Topology-based null point detection")
    print("• Magnetopause surface search algorithm")
    print("• Bisection method for precise location")
    print("• Field structure analysis around null points")
    print("• Interactive 3D visualization")
    
    print("\n⚙️ NULL POINT TYPES DETECTED:")
    print("• Topology transitions: open_north ↔ open_south")
    print("• IMF boundaries: closed → imf_draped")
    print("• Mixed transitions: various combinations")
    print("• All locations refined for minimum field strength")

