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
        GSM coordinate meshgrids from gmComputeSurface
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
    
    # Convert GSM to Geocentric Solar Magnetospheric if needed
    # (Note: GSM and Geocentric Solar Magnetospheric are the same coordinate system)
    # GSM = Geocentric Solar Magnetospheric
    
    # Convert units to meters
    if units == 'Re':
        scale_factor = earth_radius_m
    elif units == 'km':
        scale_factor = 1000
    else:  # assume meters
        scale_factor = 1
    
    x_m = x_grid * scale_factor
    y_m = y_grid * scale_factor  
    z_m = z_grid * scale_factor
    
    # Create metadata
    metadata = {
        'created': datetime.now().isoformat(),
        'coordinate_system': 'GSM (Geocentric Solar Magnetospheric)',
        'units': 'meters',
        'original_units': units,
        'grid_size': f"{x_grid.shape[0]}x{x_grid.shape[1]}",
        'earth_radius_m': earth_radius_m
    }
    
    if timestamp:
        metadata['computation_time'] = str(timestamp)
    if model_params:
        metadata['model_parameters'] = model_params
    
    # Write OBJ file
    _write_obj_file(x_m, y_m, z_m, f"{filename}.obj", metadata)
    
    # Write OpenSpace asset file
    _write_openspace_asset(filename, metadata)
    
    print(f"OpenSpace files created: {filename}.obj, {filename}.asset")


def _write_obj_file(x_grid, y_grid, z_grid, obj_filename, metadata):
    """Write OBJ mesh file with metadata"""
    
    rows, cols = x_grid.shape
    
    with open(obj_filename, 'w') as f:
        # Write metadata as comments
        f.write("# Bow Shock Surface Mesh\n")
        f.write(f"# Created: {metadata['created']}\n")
        f.write(f"# Coordinate System: {metadata['coordinate_system']}\n")
        f.write(f"# Units: {metadata['units']}\n")
        f.write(f"# Grid Size: {metadata['grid_size']}\n")
        
        if 'computation_time' in metadata:
            f.write(f"# Computation Time: {metadata['computation_time']}\n")
        
        f.write("#\n")
        
        # Write vertices
        for i in range(rows):
            for j in range(cols):
                f.write(f"v {x_grid[i,j]:.6e} {y_grid[i,j]:.6e} {z_grid[i,j]:.6e}\n")
        
        # Write faces (triangles)
        # Each grid cell creates 2 triangles
        for i in range(rows - 1):
            for j in range(cols - 1):
                # Vertex indices (OBJ uses 1-based indexing)
                v1 = i * cols + j + 1
                v2 = i * cols + (j + 1) + 1
                v3 = (i + 1) * cols + j + 1
                v4 = (i + 1) * cols + (j + 1) + 1
                
                # Two triangles per grid cell
                f.write(f"f {v1} {v2} {v3}\n")
                f.write(f"f {v2} {v4} {v3}\n")


def _write_openspace_asset(filename, metadata):
    """Write OpenSpace asset file"""
    
    asset_content = f'''-- Bow Shock Surface Asset
-- {metadata['created']}
-- {metadata['coordinate_system']}

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


# Example usage:
# x, y, z = gmComputeSurface(...)  # Your existing function
# gmSaveOpenSpace(x, y, z, "bowshock_20240101", 
#                timestamp="2024-01-01T12:00:00Z",
#                model_params={"solar_wind_pressure": 2.1, "imf_bz": -5.0})

