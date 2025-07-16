import numpy as np
from datetime import datetime, timedelta
import json

def gmSaveOpenSpace(x_grid, y_grid, z_grid, filename, 
                   timestamp=None, duration_minutes=4, model_params=None, 
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
        Start time for the bow shock surface display
    duration_minutes : int, default 4
        How long the surface should be visible (in minutes)
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
    
    # Create metadata
    metadata = {
        'created': datetime.now().isoformat(),
        'coordinate_system': 'GSM (Geocentric Solar Magnetospheric)',
        'units': 'meters',
        'original_units': units,
        'grid_size': f"{x_array.shape[0]}x{x_array.shape[1]}",
        'valid_points': int(np.sum(valid_mask)),
        'total_points': int(valid_mask.size),
        'earth_radius_m': earth_radius_m,
        'start_time': start_time,
        'end_time': end_time,
        'duration_minutes': duration_minutes
    }
    
    if model_params:
        metadata['model_parameters'] = model_params
    
    # Write OBJ file
    _write_obj_file(x_m, y_m, z_m, valid_mask, f"{filename}.obj", metadata)
    
    # Write OpenSpace asset file
    _write_openspace_asset(filename, metadata)
    
    print(f"OpenSpace files created: {filename}.obj, {filename}.asset")
    print(f"Time frame: {start_time.strftime('%Y %b %d %H:%M:%S')} to {end_time.strftime('%Y %b %d %H:%M:%S')}")


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


def _write_openspace_asset(filename, metadata):
    """Write OpenSpace asset file"""
    
    start_time_str = metadata['start_time'].strftime('%Y %b %d %H:%M:%S')
    end_time_str = metadata['end_time'].strftime('%Y %b %d %H:%M:%S')
    identifier = f"BowShock_{filename.replace('/', '_').replace('.', '_')}"
    
    asset_content = f'''-- Bow Shock Surface Asset
-- Created: {metadata['created']}
-- Coordinate System: {metadata['coordinate_system']}
-- Valid points: {metadata['valid_points']}/{metadata['total_points']}
-- Time Frame: {start_time_str} to {end_time_str}

local transforms = asset.require("scene/solarsystem/planets/earth/transforms_gsm_sm")

local bowShock = {{
    Identifier = "{identifier}",
    Parent = transforms.GeocentricSolarMagnetospheric.Identifier,
    TimeFrame = {{
        Type = "TimeFrameInterval",
        Start = "{start_time_str}",
        End = "{end_time_str}"
    }},
    Renderable = {{
        Type = "RenderableModel",
        GeometryFile = asset.resource("{filename}.obj"),
        ModelScale = 1.0, -- Already in meters
        EnableDepthTest = true,
        EnableFaceCulling = false,
        Enabled = true,
        Opacity = 0.7,
        -- Coloring options
        PerformShading = true,
        AmbientIntensity = 0.2,
        DiffuseIntensity = 0.8,
        SpecularIntensity = 0.0
    }},
    GUI = {{
        Name = "Bow Shock Surface",
        Path = "/Solar System/Earth/Magnetosphere",
        Description = "Bow shock surface computed from magnetospheric model"
    }}
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
# gmSaveOpenSpace(x, y, z, "BS2024-05-10-17-56-00", 
#                timestamp="2024-05-10T17:56:00",
#                duration_minutes=4,
#                model_params={"solar_wind_pressure": 2.1, "imf_bz": -5.0})

