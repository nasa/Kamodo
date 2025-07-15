import numpy as np
from datetime import datetime
import json

def gmSaveOpenSpace(x_grid, y_grid, z_grid, filename, 
                   timestamp=None, model_params=None, 
                   units='Re', earth_radius_m=6.371e6):
    """
    Save bow shock surface data in OpenSpace-compatible format
    
    Parameters:
    -----------
    x_grid, y_grid, z_grid : array-like (21x21)
        GSM coordinate meshgrids from gmComputeSurface (may contain None values)
    filename : str
        Output filename (without extension)
    timestamp : str or datetime, optional
        Time of computation
    model_params : dict, optional
        Model parameters used in computation
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
    
    # Create metadata
    metadata = {
        'created': datetime.now().isoformat(),
        'coordinate_system': 'GSM (Geocentric Solar Magnetospheric)',
        'units': 'meters',
        'original_units': units,
        'grid_size': f"{x_array.shape[0]}x{x_array.shape[1]}",
        'valid_points': int(np.sum(valid_mask)),
        'total_points': int(valid_mask.size),
        'earth_radius_m': earth_radius_m
    }
    
    if timestamp:
        metadata['computation_time'] = str(timestamp)
    if model_params:
        metadata['model_parameters'] = model_params
    
    # Write OBJ file
    _write_obj_file(x_m, y_m, z_m, valid_mask, f"{filename}.obj", metadata)
    
    # Write OpenSpace asset file
    _write_openspace_asset(filename, metadata)
    
    print(f"OpenSpace files created: {filename}.obj, {filename}.asset")


def _write_obj_file(x_grid, y_grid, z_grid, valid_mask, obj_filename, metadata):
    """Write OBJ mesh file with metadata, handling invalid points"""
    
    rows, cols = x_grid.shape
    
    with open(obj_filename, 'w') as f:
        # Write metadata as comments
        f.write("# Bow Shock Surface Mesh\n")
        f.write(f"# Created: {metadata['created']}\n")
        f.write(f"# Coordinate System: {metadata['coordinate_system']}\n")
        f.write(f"# Units: {metadata['units']}\n")
        f.write(f"# Grid Size: {metadata['grid_size']}\n")
        f.write(f"# Valid Points: {metadata['valid_points']}/{metadata['total_points']}\n")
        
        if 'computation_time' in metadata:
            f.write(f"# Computation Time: {metadata['computation_time']}\n")
        
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


def _write_openspace_asset(filename, metadata):
    """Write OpenSpace asset file"""
    
    asset_content = f'''-- Bow Shock Surface Asset
-- {metadata['created']}
-- {metadata['coordinate_system']}
-- Valid points: {metadata['valid_points']}/{metadata['total_points']}

local bowShock = {{
    Identifier = "BowShock_{filename.replace('/', '_').replace('.', '_')}",
    Renderable = {{
        Type = "RenderableModel",
        GeometryFile = "{filename}.obj",
        ModelScale = 1.0, -- Already in meters
        EnableDepthTest = true,
        Enabled = true,
        Opacity = 0.7,
        -- Coloring options
        PerformShading = true,
        AmbientIntensity = 0.2,
        DiffuseIntensity = 0.8,
        SpecularIntensity = 0.0
    }},
    Transform = {{
        Translation = {{
            Type = "StaticTranslation",
            Position = {{ 0, 0, 0 }} -- Relative to Earth center
        }}
    }},
    GUI = {{
        Name = "Bow Shock Surface",
        Path = "/Solar System/Earth/Magnetosphere",
        Description = "Bow shock surface computed from magnetospheric model"
    }},
    -- Metadata
    Tag = {{ "earth_magnetosphere", "bow_shock" }}
}}

-- Asset metadata
bowShock.Meta = {{
    Name = "Bow Shock Surface",
    Version = "1.0",
    Description = "Magnetospheric bow shock surface",
    CreatedBy = "Kamodo gmSaveOpenSpace function",'''

    if 'model_parameters' in metadata:
        asset_content += f'''
    ModelParameters = {json.dumps(metadata['model_parameters'], indent=8).replace('{', '{{').replace('}', '}}')}'''

    asset_content += f'''
}}

asset.onInitialize(function()
    openspace.addSceneGraphNode(bowShock)
end)

asset.onDeinitialize(function()
    openspace.removeSceneGraphNode(bowShock)
end)

asset.export(bowShock)
'''

    with open(f"{filename}.asset", 'w') as f:
        f.write(asset_content)

