import numpy as np
from datetime import datetime, timedelta
import json
import os
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def gmSaveOpenSpace(x_grid, y_grid, z_grid, filename, surface_type,
                   timestamp=None, duration_minutes=4, model_params=None,
                   color=None, color_values=None, colormap='plasma',
                   color_min=None, color_max=None, 
                   units='Re', earth_radius_m=6.371e6):
    """
    Save 3D surface data in OpenSpace-compatible format
    
    Parameters:
    -----------
    x_grid, y_grid, z_grid : array-like (NxM)
        GSM coordinate meshgrids from surface computation (may contain None values)
    filename : str
        Output filename (without extension)
    surface_type : str
        Type of surface (e.g., 'bowshock', 'magnetopause', 'slice', 'boundary')
    timestamp : str or datetime, optional
        Start time for the surface display
    duration_minutes : int, default 4
        How long the surface should be visible (in minutes)
    model_params : dict, optional
        Model parameters used in computation
    color : tuple or list, optional
        RGB color for the entire surface, e.g., (1.0, 0.5, 0.0) for orange
        If provided, overrides color_values
    color_values : array-like, same shape as x_grid, optional
        Scalar values to map to colors (like pressure, temperature, density)
    colormap : str, default 'plasma'
        Matplotlib colormap name to use for scalar values
    color_min, color_max : float, optional
        Min/max values for color scaling. If None, auto-determined from data
    units : str, default 'Re'
        Input units ('Re', 'km', 'm')
    earth_radius_m : float, default 6.371e6
        Earth radius in meters for unit conversion
    """
    
    # Handle None values by converting to numpy arrays with masked arrays or NaN
    x_array = np.array(x_grid, dtype=float)
    y_array = np.array(y_grid, dtype=float)
    z_array = np.array(z_grid, dtype=float)
    
    # Create mask for valid (non-None, non-NaN) values
    valid_mask = (~np.isnan(x_array)) & (~np.isnan(y_array)) & (~np.isnan(z_array))
    
    # Also check color_values validity if provided
    use_vertex_colors = False
    if color_values is not None and color is None:
        color_array = np.array(color_values, dtype=float)
        valid_mask = valid_mask & (~np.isnan(color_array))
        use_vertex_colors = True
        
        # Set up colormap
        cmap = plt.get_cmap(colormap)
        vmin = color_min if color_min is not None else np.nanmin(color_array)
        vmax = color_max if color_max is not None else np.nanmax(color_array)
        norm = Normalize(vmin=vmin, vmax=vmax)
        
        # Generate RGB colors for each vertex
        vertex_colors = {}
        for i in range(x_array.shape[0]):
            for j in range(x_array.shape[1]):
                if valid_mask[i, j]:
                    # Get normalized value and convert to RGB
                    value = color_array[i, j]
                    rgb = cmap(norm(value))[:3]  # Exclude alpha
                    vertex_colors[(i, j)] = rgb
                    
        color_info = {
            'type': 'vertex',
            'colormap': colormap,
            'min': float(vmin),
            'max': float(vmax),
            'data_type': 'Custom data'
        }
    else:
        # Uniform color for the whole surface
        color_info = {
            'type': 'uniform',
            'color': color
        }
    
    if not np.any(valid_mask):
        raise ValueError("No valid (non-None) data points found in grids")
    
    print(f"Found {np.sum(valid_mask)} valid points out of {valid_mask.size} total points")
    
    # Convert units to meters
    if units == 'Re':
        scale_factor = earth_radius_m
    elif units == 'km':
        scale_factor = 1000
    else:  # assume meters
        scale_factor = 1
    
    x_m = x_array * scale_factor
    y_m = y_array * scale_factor  
    z_m = z_array * scale_factor
    
    # Handle timestamp and duration
    if timestamp is None:
        start_time = datetime.now()
    elif isinstance(timestamp, str):
        # Try to parse common formats
        try:
            start_time = datetime.fromisoformat(timestamp.replace('Z', '+00:00'))
        except:
            # If parsing fails, use current time
            start_time = datetime.now()
            print(f"Warning: Could not parse timestamp '{timestamp}', using current time")
    else:
        start_time = timestamp
    
    end_time = start_time + timedelta(minutes=duration_minutes)
    
    # Create metadata with surface type information
    metadata = {
        'created': datetime.now().isoformat(),
        'surface_type': surface_type.lower(),
        'coordinate_system': 'GSM (Geocentric Solar Magnetospheric)',
        'units': 'meters',
        'original_units': units,
        'grid_size': f"{x_array.shape[0]}x{x_array.shape[1]}",
        'valid_points': int(np.sum(valid_mask)),
        'total_points': int(valid_mask.size),
        'earth_radius_m': earth_radius_m,
        'start_time': start_time,
        'end_time': end_time,
        'duration_minutes': duration_minutes,
        'color_info': color_info
    }
    
    if model_params:
        metadata['model_parameters'] = model_params
    
    # Write OBJ file
    if use_vertex_colors:
        _write_colored_obj_file(x_m, y_m, z_m, valid_mask, vertex_colors, f"{filename}.obj", metadata)
        _write_mtl_file(filename, metadata)
    else:
        _write_obj_file(x_m, y_m, z_m, valid_mask, f"{filename}.obj", metadata)
    
    # Write OpenSpace asset file
    _write_openspace_asset(filename, metadata)
    
    # Generate a colormap legend if using vertex colors
    if use_vertex_colors:
        _save_colormap_legend(filename, cmap, vmin, vmax, color_info.get('data_type', 'Value'))
    
    print(f"OpenSpace files created: {filename}.obj, {filename}.asset")
    if use_vertex_colors:
        print(f"Also created: {filename}.mtl (material file) and {filename}_colormap.png (legend)")
    print(f"Time frame: {start_time.strftime('%Y %b %d %H:%M:%S')} to {end_time.strftime('%Y %b %d %H:%M:%S')}")


def gmSaveOpenSpaceIndexed(x_array, y_array, z_array, i_indices, j_indices, k_indices, 
                          filename, surface_type, timestamp=None, duration_minutes=4, 
                          model_params=None, color=None, color_values=None, 
                          colormap='plasma', color_min=None, color_max=None, 
                          units='Re', earth_radius_m=6.371e6):
    """
    Save 3D surface data with indexed connectivity in OpenSpace-compatible format
    
    Parameters:
    -----------
    x_array, y_array, z_array : array-like (N,)
        1D arrays of GSM coordinate positions (may contain None values)
    i_indices, j_indices, k_indices : array-like (M,)
        1D integer arrays defining triangle vertices (0-based indices into coordinate arrays)
        Each element represents one vertex of a triangle, so M triangles total
    filename : str
        Output filename (without extension)
    surface_type : str
        Type of surface (e.g., 'bowshock', 'magnetopause', 'slice', 'boundary')
    timestamp : str or datetime, optional
        Start time for the surface display
    duration_minutes : int, default 4
        How long the surface should be visible (in minutes)
    model_params : dict, optional
        Model parameters used in computation
    color : tuple or list, optional
        RGB color for the entire surface, e.g., (1.0, 0.5, 0.0) for orange
        If provided, overrides color_values
    color_values : array-like (N,), optional
        Scalar values to map to colors (same length as coordinate arrays)
    colormap : str, default 'plasma'
        Matplotlib colormap name to use for scalar values
    color_min, color_max : float, optional
        Min/max values for color scaling. If None, auto-determined from data
    units : str, default 'Re'
        Input units ('Re', 'km', 'm')
    earth_radius_m : float, default 6.371e6
        Earth radius in meters for unit conversion
    """
    
    # Convert to numpy arrays
    x_coord = np.array(x_array, dtype=float)
    y_coord = np.array(y_array, dtype=float)
    z_coord = np.array(z_array, dtype=float)
    
    # Convert indices to numpy arrays
    i_idx = np.array(i_indices, dtype=int)
    j_idx = np.array(j_indices, dtype=int)
    k_idx = np.array(k_indices, dtype=int)
    
    # Validate that coordinate arrays are same length
    if not (len(x_coord) == len(y_coord) == len(z_coord)):
        raise ValueError("Coordinate arrays must be the same length")
    
    # Validate that index arrays are same length (each represents triangles)
    if not (len(i_idx) == len(j_idx) == len(k_idx)):
        raise ValueError("Index arrays must be the same length")
    
    num_triangles = len(i_idx)  # Each element is one triangle
    
    # Create mask for valid coordinates
    valid_coords = (~np.isnan(x_coord)) & (~np.isnan(y_coord)) & (~np.isnan(z_coord))
    
    # Handle color values
    use_vertex_colors = False
    if color_values is not None and color is None:
        color_array = np.array(color_values, dtype=float)
        if len(color_array) != len(x_coord):
            raise ValueError("color_values must be same length as coordinate arrays")
        
        # Update valid mask to include color validity
        valid_coords = valid_coords & (~np.isnan(color_array))
        use_vertex_colors = True
        
        # Set up colormap
        cmap = plt.get_cmap(colormap)
        vmin = color_min if color_min is not None else np.nanmin(color_array[valid_coords])
        vmax = color_max if color_max is not None else np.nanmax(color_array[valid_coords])
        norm = Normalize(vmin=vmin, vmax=vmax)
        
        color_info = {
            'type': 'vertex',
            'colormap': colormap,
            'min': float(vmin),
            'max': float(vmax),
            'data_type': 'Custom data'
        }
    else:
        color_info = {
            'type': 'uniform',
            'color': color
        }
    
    if not np.any(valid_coords):
        raise ValueError("No valid coordinate data found")
    
    # Filter out triangles that reference invalid vertices
    valid_triangles = []
    
    for t in range(num_triangles):
        v1 = i_idx[t]  # First vertex of triangle t
        v2 = j_idx[t]  # Second vertex of triangle t  
        v3 = k_idx[t]  # Third vertex of triangle t
        
        # Check if all vertex indices are valid and within bounds
        if (v1 >= 0 and v2 >= 0 and v3 >= 0 and
            v1 < len(x_coord) and v2 < len(x_coord) and v3 < len(x_coord) and
            valid_coords[v1] and valid_coords[v2] and valid_coords[v3]):
            valid_triangles.append([v1, v2, v3])
    
    if len(valid_triangles) == 0:
        raise ValueError("No valid triangles found")
    
    valid_triangles = np.array(valid_triangles)
    
    print(f"Found {np.sum(valid_coords)} valid vertices out of {len(x_coord)} total")
    print(f"Found {len(valid_triangles)} valid triangles out of {num_triangles} total")
    
    # Convert units to meters
    if units == 'Re':
        scale_factor = earth_radius_m
    elif units == 'km':
        scale_factor = 1000
    else:  # assume meters
        scale_factor = 1
    
    x_m = x_coord * scale_factor
    y_m = y_coord * scale_factor  
    z_m = z_coord * scale_factor
    
    # Handle timestamp and duration
    if timestamp is None:
        start_time = datetime.now()
    elif isinstance(timestamp, str):
        try:
            start_time = datetime.fromisoformat(timestamp.replace('Z', '+00:00'))
        except:
            start_time = datetime.now()
            print(f"Warning: Could not parse timestamp '{timestamp}', using current time")
    else:
        start_time = timestamp
    
    end_time = start_time + timedelta(minutes=duration_minutes)
    
    # Create metadata
    metadata = {
        'created': datetime.now().isoformat(),
        'surface_type': surface_type.lower(),
        'coordinate_system': 'GSM (Geocentric Solar Magnetospheric)',
        'units': 'meters',
        'original_units': units,
        'data_structure': 'indexed',
        'num_vertices': len(x_coord),
        'num_triangles': num_triangles,
        'valid_vertices': int(np.sum(valid_coords)),
        'valid_triangles': len(valid_triangles),
        'earth_radius_m': earth_radius_m,
        'start_time': start_time,
        'end_time': end_time,
        'duration_minutes': duration_minutes,
        'color_info': color_info
    }
    
    if model_params:
        metadata['model_parameters'] = model_params
    
    # Write files
    if use_vertex_colors:
        _write_indexed_colored_obj_file(x_m, y_m, z_m, valid_triangles, valid_coords, 
                                       color_array, norm, cmap, f"{filename}.obj", metadata)
        _write_mtl_file(filename, metadata)
        _save_colormap_legend(filename, cmap, vmin, vmax, color_info.get('data_type', 'Value'))
    else:
        _write_indexed_obj_file(x_m, y_m, z_m, valid_triangles, valid_coords, 
                               f"{filename}.obj", metadata)
    
    _write_openspace_asset(filename, metadata)
    
    print(f"OpenSpace files created: {filename}.obj, {filename}.asset")
    if use_vertex_colors:
        print(f"Also created: {filename}.mtl and {filename}_colormap.png")
    print(f"Time frame: {start_time.strftime('%Y %b %d %H:%M:%S')} to {end_time.strftime('%Y %b %d %H:%M:%S')}")


def _get_surface_display_info(surface_type):
    """Get display information based on surface type"""
    
    surface_configs = {
        'bowshock': {
            'name': 'Bow Shock Surface',
            'description': 'Magnetospheric bow shock surface',
            'gui_path': '/Solar System/Earth/Magnetosphere',
            'tags': ['earth_magnetosphere', 'bow_shock']
        },
        'magnetopause': {
            'name': 'Magnetopause Surface',
            'description': 'Magnetospheric boundary surface',
            'gui_path': '/Solar System/Earth/Magnetosphere',
            'tags': ['earth_magnetosphere', 'magnetopause']
        },
        'slice': {
            'name': '3D Data Slice',
            'description': 'Cross-section through 3D data',
            'gui_path': '/Solar System/Earth/Data',
            'tags': ['earth_data', 'slice']
        },
        'boundary': {
            'name': 'Boundary Surface',
            'description': 'Generic boundary surface',
            'gui_path': '/Solar System/Earth/Boundaries',
            'tags': ['earth_boundaries', 'surface']
        },
        'plasmapause': {
            'name': 'Plasmapause Surface',
            'description': 'Inner magnetospheric plasma boundary',
            'gui_path': '/Solar System/Earth/Magnetosphere',
            'tags': ['earth_magnetosphere', 'plasmapause']
        },
        'shock': {
            'name': 'Shock Surface',
            'description': 'Plasma shock boundary',
            'gui_path': '/Solar System/Earth/Plasma',
            'tags': ['earth_plasma', 'shock']
        }
    }
    
    # Default configuration if surface type not found
    default_config = {
        'name': f'{surface_type.title()} Surface',
        'description': f'{surface_type.title()} surface from magnetospheric model',
        'gui_path': '/Solar System/Earth/Surfaces',
        'tags': ['earth_surfaces', surface_type.lower()]
    }
    
    return surface_configs.get(surface_type.lower(), default_config)


def _write_obj_file(x_grid, y_grid, z_grid, valid_mask, obj_filename, metadata):
    """Write OBJ mesh file with metadata, handling invalid points"""
    
    rows, cols = x_grid.shape
    surface_type = metadata.get('surface_type', 'surface')
    
    with open(obj_filename, 'w') as f:
        # Write metadata as comments
        f.write(f"# {surface_type.title()} Surface Mesh\n")
        f.write(f"# Created: {metadata['created']}\n")
        f.write(f"# Surface Type: {surface_type}\n")
        f.write(f"# Coordinate System: {metadata['coordinate_system']}\n")
        f.write(f"# Units: {metadata['units']}\n")
        f.write(f"# Grid Size: {metadata['grid_size']}\n")
        f.write(f"# Valid Points: {metadata['valid_points']}/{metadata['total_points']}\n")
        f.write(f"# Time Frame: {metadata['start_time']} to {metadata['end_time']}\n")
        
        if 'model_parameters' in metadata:
            f.write(f"# Model Parameters: {metadata['model_parameters']}\n")
        
        f.write("#\n")
        
        # Create vertex index mapping (only for valid points)
        vertex_map = {}
        vertex_count = 0
        
        # Write vertices (only valid ones)
        for i in range(rows):
            for j in range(cols):
                if valid_mask[i, j]:
                    vertex_count += 1
                    vertex_map[(i, j)] = vertex_count
                    f.write(f"v {x_grid[i,j]:.6e} {y_grid[i,j]:.6e} {z_grid[i,j]:.6e}\n")
        
        # Write faces (triangles) - only for cells where all 4 corners are valid
        for i in range(rows - 1):
            for j in range(cols - 1):
                # Check if all 4 corners of this grid cell are valid
                corners_valid = (valid_mask[i, j] and valid_mask[i, j+1] and 
                               valid_mask[i+1, j] and valid_mask[i+1, j+1])
                
                if corners_valid:
                    # Get vertex indices (OBJ uses 1-based indexing)
                    v1 = vertex_map[(i, j)]
                    v2 = vertex_map[(i, j+1)]
                    v3 = vertex_map[(i+1, j)]
                    v4 = vertex_map[(i+1, j+1)]
                    
                    # Two triangles per grid cell
                    f.write(f"f {v1} {v2} {v3}\n")
                    f.write(f"f {v2} {v4} {v3}\n")


def _write_colored_obj_file(x_grid, y_grid, z_grid, valid_mask, vertex_colors, obj_filename, metadata):
    """Write OBJ mesh file with vertex colors"""
    
    rows, cols = x_grid.shape
    surface_type = metadata.get('surface_type', 'surface')
    
    with open(obj_filename, 'w') as f:
        # Write metadata as comments
        f.write(f"# {surface_type.title()} Surface Mesh with Vertex Colors\n")
        f.write(f"# Created: {metadata['created']}\n")
        f.write(f"# Surface Type: {surface_type}\n")
        f.write(f"# Coordinate System: {metadata['coordinate_system']}\n")
        f.write(f"# Units: {metadata['units']}\n")
        f.write(f"# Grid Size: {metadata['grid_size']}\n")
        f.write(f"# Valid Points: {metadata['valid_points']}/{metadata['total_points']}\n")
        f.write(f"# Time Frame: {metadata['start_time']} to {metadata['end_time']}\n")
        
        if 'model_parameters' in metadata:
            f.write(f"# Model Parameters: {metadata['model_parameters']}\n")
        
        # Reference material library
        mtl_filename = obj_filename.replace('.obj', '.mtl')
        f.write(f"mtllib {os.path.basename(mtl_filename)}\n\n")
        
        # Create vertex index mapping (only for valid points)
        vertex_map = {}
        vertex_count = 0
        
        # Write vertices with vertex colors
        for i in range(rows):
            for j in range(cols):
                if valid_mask[i, j]:
                    vertex_count += 1
                    vertex_map[(i, j)] = vertex_count
                    
                    # Get RGB color for this vertex
                    r, g, b = vertex_colors.get((i, j), (1.0, 1.0, 1.0))
                    
                    # Write vertex position and color
                    f.write(f"v {x_grid[i,j]:.6e} {y_grid[i,j]:.6e} {z_grid[i,j]:.6e} {r:.6f} {g:.6f} {b:.6f}\n")
        
        # Set material
        f.write(f"\nusemtl {surface_type.title()}Material\n")
        
        # Write faces (triangles) - only for cells where all 4 corners are valid
        f.write("\n# Faces\n")
        for i in range(rows - 1):
            for j in range(cols - 1):
                # Check if all 4 corners of this grid cell are valid
                corners_valid = (valid_mask[i, j] and valid_mask[i, j+1] and 
                               valid_mask[i+1, j] and valid_mask[i+1, j+1])
                
                if corners_valid:
                    # Get vertex indices (OBJ uses 1-based indexing)
                    v1 = vertex_map[(i, j)]
                    v2 = vertex_map[(i, j+1)]
                    v3 = vertex_map[(i+1, j)]
                    v4 = vertex_map[(i+1, j+1)]
                    
                    # Two triangles per grid cell
                    f.write(f"f {v1} {v2} {v3}\n")
                    f.write(f"f {v2} {v4} {v3}\n")


def _write_indexed_obj_file(x_coords, y_coords, z_coords, triangles, valid_coords, obj_filename, metadata):
    """Write OBJ file for indexed vertex data"""
    
    surface_type = metadata.get('surface_type', 'surface')
    
    with open(obj_filename, 'w') as f:
        # Write metadata as comments
        f.write(f"# {surface_type.title()} Surface Mesh (Indexed)\n")
        f.write(f"# Created: {metadata['created']}\n")
        f.write(f"# Surface Type: {surface_type}\n")
        f.write(f"# Coordinate System: {metadata['coordinate_system']}\n")
        f.write(f"# Units: {metadata['units']}\n")
        f.write(f"# Data Structure: {metadata['data_structure']}\n")
        f.write(f"# Vertices: {metadata['valid_vertices']}/{metadata['num_vertices']}\n")
        f.write(f"# Triangles: {metadata['valid_triangles']}/{metadata['num_triangles']}\n")
        f.write(f"# Time Frame: {metadata['start_time']} to {metadata['end_time']}\n")
        
        if 'model_parameters' in metadata:
            f.write(f"# Model Parameters: {metadata['model_parameters']}\n")
        
        f.write("#\n")
        
        # Create mapping from original indices to OBJ indices (1-based)
        vertex_map = {}
        obj_vertex_count = 0
        
        # Write vertices (only valid ones)
        for i, (valid, x, y, z) in enumerate(zip(valid_coords, x_coords, y_coords, z_coords)):
            if valid:
                obj_vertex_count += 1
                vertex_map[i] = obj_vertex_count
                f.write(f"v {x:.6e} {y:.6e} {z:.6e}\n")
        
        # Write faces using the vertex mapping
        f.write("\n# Faces\n")
        for triangle in triangles:
            v1, v2, v3 = triangle
            if all(idx in vertex_map for idx in triangle):
                obj_v1 = vertex_map[v1]
                obj_v2 = vertex_map[v2]
                obj_v3 = vertex_map[v3]
                f.write(f"f {obj_v1} {obj_v2} {obj_v3}\n")


def _write_indexed_colored_obj_file(x_coords, y_coords, z_coords, triangles, valid_coords, 
                                   color_values, norm, cmap, obj_filename, metadata):
    """Write OBJ file for indexed vertex data with colors"""
    
    surface_type = metadata.get('surface_type', 'surface')
    
    with open(obj_filename, 'w') as f:
        # Write metadata as comments
        f.write(f"# {surface_type.title()} Surface Mesh (Indexed, Colored)\n")
        f.write(f"# Created: {metadata['created']}\n")
        f.write(f"# Surface Type: {surface_type}\n")
        f.write(f"# Coordinate System: {metadata['coordinate_system']}\n")
        f.write(f"# Units: {metadata['units']}\n")
        f.write(f"# Data Structure: {metadata['data_structure']}\n")
        f.write(f"# Vertices: {metadata['valid_vertices']}/{metadata['num_vertices']}\n")
        f.write(f"# Triangles: {metadata['valid_triangles']}/{metadata['num_triangles']}\n")
        f.write(f"# Time Frame: {metadata['start_time']} to {metadata['end_time']}\n")
        
        if 'model_parameters' in metadata:
            f.write(f"# Model Parameters: {metadata['model_parameters']}\n")
        
        # Reference material library
        mtl_filename = obj_filename.replace('.obj', '.mtl')
        f.write(f"mtllib {os.path.basename(mtl_filename)}\n\n")
        
        # Create mapping from original indices to OBJ indices (1-based)
        vertex_map = {}
        obj_vertex_count = 0
        
        # Write vertices with colors (only valid ones)
        for i, (valid, x, y, z) in enumerate(zip(valid_coords, x_coords, y_coords, z_coords)):
            if valid:
                obj_vertex_count += 1
                vertex_map[i] = obj_vertex_count
                
                # Get color for this vertex
                color_val = color_values[i]
                rgb = cmap(norm(color_val))[:3]  # Exclude alpha
                
                f.write(f"v {x:.6e} {y:.6e} {z:.6e} {rgb[0]:.6f} {rgb[1]:.6f} {rgb[2]:.6f}\n")
        
        # Set material
        f.write(f"\nusemtl {surface_type.title()}Material\n")
        
        # Write faces using the vertex mapping
        f.write("\n# Faces\n")
        for triangle in triangles:
            v1, v2, v3 = triangle
            if all(idx in vertex_map for idx in triangle):
                obj_v1 = vertex_map[v1]
                obj_v2 = vertex_map[v2]
                obj_v3 = vertex_map[v3]
                f.write(f"f {obj_v1} {obj_v2} {obj_v3}\n")


def _write_mtl_file(filename, metadata):
    """Write material file for colored OBJ"""
    
    surface_type = metadata.get('surface_type', 'surface')
    mtl_filename = f"{filename}.mtl"
    
    with open(mtl_filename, 'w') as f:
        f.write(f"# {surface_type.title()} Material\n")
        f.write(f"# Created: {metadata['created']}\n\n")
        
        f.write(f"newmtl {surface_type.title()}Material\n")
        f.write("Ka 0.2 0.2 0.2\n")  # Ambient color
        f.write("Kd 0.8 0.8 0.8\n")  # Diffuse color
        f.write("Ks 0.0 0.0 0.0\n")  # Specular color
        f.write("d 0.7\n")           # Transparency (0.7 opacity)
        f.write("illum 2\n")         # Illumination model


def _save_colormap_legend(filename, cmap, vmin, vmax, title='Value'):
    """Generate a colormap legend image"""
    
    fig, ax = plt.subplots(figsize=(6, 1))
    
    # Create a colorbar
    norm = Normalize(vmin=vmin, vmax=vmax)
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), 
                      cax=ax, orientation='horizontal')
    
    cb.set_label(title)
    plt.savefig(f"{filename}_colormap.png", dpi=150, bbox_inches='tight')
    plt.close()


def _write_openspace_asset(filename, metadata):
    """Write OpenSpace asset file (updated to handle both grid and indexed data)"""
    
    start_time_str = metadata['start_time'].strftime('%Y %b %d %H:%M:%S')
    end_time_str = metadata['end_time'].strftime('%Y %b %d %H:%M:%S')
    surface_type = metadata.get('surface_type', 'surface')
    
    # Get surface-specific display information
    display_info = _get_surface_display_info(surface_type)
    
    identifier = f"{surface_type.title()}_{filename.replace('/', '_').replace('.', '_')}"
    
    # Build color specification
    color_info = metadata.get('color_info', {})
    color_spec = ""
    
    if color_info.get('type') == 'uniform' and color_info.get('color'):
        r, g, b = color_info['color']
        color_spec = f'''
        Color = {{ {r:.3f}, {g:.3f}, {b:.3f} }},'''
    
    # Handle vertex colors by enabling them in OpenSpace
    if color_info.get('type') == 'vertex':
        color_spec = '''
        UseVertexColors = true,'''
    
    # Build point info string based on data structure
    data_structure = metadata.get('data_structure', 'grid')
    if data_structure == 'indexed':
        point_info = f"Vertices: {metadata['valid_vertices']}/{metadata['num_vertices']}, Triangles: {metadata['valid_triangles']}/{metadata['num_triangles']}"
    else:
        point_info = f"Valid points: {metadata['valid_points']}/{metadata['total_points']}"
    
    asset_content = f'''-- {display_info['name']} Asset
-- Created: {metadata['created']}
-- Surface Type: {surface_type}
-- Data Structure: {data_structure}
-- Coordinate System: {metadata['coordinate_system']}
-- {point_info}
-- Time Frame: {start_time_str} to {end_time_str}'''

    if color_info.get('type') == 'vertex':
        asset_content += f'''
-- Colored by: {color_info.get('data_type', 'custom data')}
-- Colormap: {color_info.get('colormap')}
-- Range: [{color_info.get('min')}, {color_info.get('max')}]'''

    asset_content += f'''

local transforms = asset.require("scene/solarsystem/planets/earth/transforms_gsm_sm")

local {surface_type}Surface = {{
    Identifier = "{identifier}",
    Parent = transforms.GeocentricSolarMagnetospheric.Identifier,
    TimeFrame = {{
        Type = "TimeFrameInterval",
        Start = "{start_time_str}",
        End = "{end_time_str}"
    }},
    Renderable = {{
        Type = "RenderableModel",
        GeometryFile = asset.resource("{filename}.obj"),'''
    
    # Add material file if vertex colors are used
    if color_info.get('type') == 'vertex':
        asset_content += f'''
        MaterialFile = asset.resource("{filename}.mtl"),'''
    
    asset_content += f'''
        ModelScale = 1.0, -- Already in meters
        EnableDepthTest = true,
        EnableFaceCulling = false,
        Enabled = true,
        Opacity = 0.7,{color_spec}
        -- Coloring options
        PerformShading = true,
        AmbientIntensity = 0.2,
        DiffuseIntensity = 0.8,
        SpecularIntensity = 0.0
    }},
    GUI = {{
        Name = "{display_info['name']}",
        Path = "{display_info['gui_path']}",
        Description = "{display_info['description']}"
    }},
    -- Metadata
    Tag = {display_info['tags']}
}}

asset.onInitialize(function()
    openspace.addSceneGraphNode({surface_type}Surface)
end)

asset.onDeinitialize(function()
    openspace.removeSceneGraphNode({surface_type}Surface)
end)

asset.export({surface_type}Surface)
'''

    with open(f"{filename}.asset", 'w') as f:
        f.write(asset_content)

