Vector Field Tracer Documentation
Overview
The vectorfieldtracer.py module provides high-performance 3D vector field line tracing capabilities for Kamodo magnetospheric models. It supports tracing magnetic field lines, velocity streamlines, current density pathways, and electric field lines with adaptive step sizing and intelligent stopping criteria.

All coordinates are in Earth radii (RE) units.

Installation and Import
python
Copy code
import kamodo_ccmc.tools.vectorfieldtracer as vft
Key Features
✅ Multi-field support: Magnetic, velocity, current, and electric fields
✅ Adaptive step sizing: Automatic adjustment based on field strength
✅ Intelligent stopping criteria: Field-type aware boundary detection
✅ Bidirectional tracing: Forward, backward, or combined tracing
✅ Automatic bounds detection: Uses Kamodo model boundaries
✅ High performance: Optimized numerical integration
✅ Comprehensive diagnostics: Built-in testing and inspection tools
Quick Start
python
Copy code
import kamodo_ccmc.tools.vectorfieldtracer as vft

# Load your Kamodo object
ko = your_kamodo_model  # Your loaded Kamodo object

# Create a magnetic field tracer
tracer = vft.create_magnetic_tracer(ko)

# Define starting points (in Earth radii)
start_points_re = [
    (3.0, 0.0, 0.0),    # Dayside equatorial
    (-10.0, 0.0, 0.0),  # Nightside equatorial
    (0.0, 5.0, 0.0),    # Dawn flank
    (0.0, 0.0, 5.0)     # North pole
]

# Trace field lines
time_hours = 12.5  # Hours since midnight of first day
traces = tracer.trace_multiple_lines(start_points_re, time_hours)

# Get statistics
stats = tracer.get_field_statistics(traces)
print(f"Successfully traced: {stats['n_successful_traces']} field lines")
print(f"Mean length: {stats['mean_length_re']:.2f} RE")
Core Classes and Functions
KamodoVectorFieldTracer Class
The main class for tracing vector field lines through Kamodo models.

Constructor
python
Copy code
tracer = vft.KamodoVectorFieldTracer(
    kamodo_object=ko,                     # Your Kamodo model
    vector_components=('B_x', 'B_y', 'B_z'), # Vector field components
    coord_system='GSM',                   # Coordinate system
    field_type='magnetic',                # Field type
    auto_bounds=True,                     # Auto-detect model bounds
    verbose=False                         # Print detailed info
)
Parameters
Parameter	Type	Default	Description
kamodo_object	Kamodo	Required	Kamodo object containing vector field data
vector_components	tuple	Required	Names of X, Y, Z components (e.g., ('B_x', 'B_y', 'B_z'))
coord_system	str	'GSM'	Coordinate system ('GSM', 'GSE', 'SM', etc.)
field_type	str	'magnetic'	Field type: 'magnetic', 'velocity', 'current', 'electric'
auto_bounds	bool	True	Whether to automatically detect model boundaries
verbose	bool	False	Whether to print detailed information
Key Methods
trace_vector_line(start_point_re, time_hours, **kwargs)
Trace a single vector field line from a starting point.

python
Copy code
trace_result = tracer.trace_vector_line(
    start_point_re=(3.0, 0.0, 0.0),  # Starting point in RE
    time_hours=12.0,                 # Time since midnight (hours)
    max_steps=10000,                 # Maximum integration steps
    step_size_re=0.1,               # Step size in RE
    adaptive_step=True,              # Use adaptive stepping
    direction='both'                 # 'forward', 'backward', or 'both'
)
Returns: Dictionary with trace results containing:

'forward': Forward direction trace
'backward': Backward direction trace
'combined': Combined bidirectional trace (if direction='both')
trace_multiple_lines(start_points_re, time_hours, **kwargs)
Trace multiple field lines efficiently.

python
Copy code
traces = tracer.trace_multiple_lines(
    start_points_re=[(3.0, 0.0, 0.0), (0.0, 3.0, 0.0)],  # List of starting points
    time_hours=12.0,
    max_steps=5000,
    direction='both'
)
Returns: List of trace result dictionaries.

get_field_statistics(trace_results)
Compute statistics for a set of traced field lines.

python
Copy code
stats = tracer.get_field_statistics(traces)
Returns: Dictionary with statistics including:

'n_successful_traces': Number of successful traces
'mean_length_re': Mean field line length in RE
'stop_reasons': Dictionary of stopping reasons and counts
Convenience Functions
Pre-configured tracer creation functions for common field types:

python
Copy code
# Magnetic field tracer (B_x, B_y, B_z)
mag_tracer = vft.create_magnetic_tracer(ko)

# Velocity field tracer (v_x, v_y, v_z)
vel_tracer = vft.create_velocity_tracer(ko)

# Current density tracer (J_x, J_y, J_z)
cur_tracer = vft.create_current_tracer(ko)

# Electric field tracer (E_x, E_y, E_z)
elec_tracer = vft.create_electric_tracer(ko)
Diagnostic Functions
inspect_kamodo_object(ko)
Inspect the structure and available variables in a Kamodo object.

python
Copy code
available_vars = vft.inspect_kamodo_object(ko)
test_kamodo_integration(ko, test_time=12.0)
Test integration with your Kamodo object and diagnose any issues.

python
Copy code
vft.test_kamodo_integration(ko, test_time=12.0)
get_model_bounds(ko, vector_components)
Get coordinate bounds from your Kamodo model.

python
Copy code
bounds = vft.get_model_bounds(ko, ('B_x', 'B_y', 'B_z'))
print(f"X range: {bounds['x_range']} RE")
print(f"Y range: {bounds['y_range']} RE")
print(f"Z range: {bounds['z_range']} RE")
Field Types and Configuration
The tracer automatically configures parameters based on field type:

Magnetic Fields (field_type='magnetic')
Integration method: Field line tracing
Default step size: 0.1 RE
Stopping criteria: Earth surface, magnetopause, weak field
Max distance: 100 RE
Velocity Fields (field_type='velocity')
Integration method: Streamline tracing
Default step size: 0.2 RE
Stopping criteria: Flow boundaries, Earth surface
Max distance: 50 RE
Current Fields (field_type='current')
Integration method: Current streamline tracing
Default step size: 0.15 RE
Stopping criteria: Current sheet boundaries
Max distance: 30 RE
Electric Fields (field_type='electric')
Integration method: Field line tracing
Default step size: 0.1 RE
Stopping criteria: Ionospheric boundaries
Max distance: 20 RE
Usage Examples
Example 1: Basic Field Line Tracing
python
Copy code
import kamodo_ccmc.tools.vectorfieldtracer as vft
import numpy as np

# Load your Kamodo object
ko = your_kamodo_model

# Create magnetic field tracer
tracer = vft.create_magnetic_tracer(ko)

# Define starting points in the equatorial plane
n_points = 8
angles = np.linspace(0, 2*np.pi, n_points, endpoint=False)
radius = 4.0  # 4 RE from Earth center

start_points_re = []
for angle in angles:
    x = radius * np.cos(angle)
    y = radius * np.sin(angle)
    z = 0.0
    start_points_re.append((x, y, z))

print(f"Tracing {len(start_points_re)} field lines from equatorial plane")

# Trace field lines
time_hours = 12.0
traces = tracer.trace_multiple_lines(
    start_points_re, 
    time_hours,
    max_steps=5000,
    step_size_re=0.1,
    direction='both',
    adaptive_step=True
)

# Analyze results
stats = tracer.get_field_statistics(traces)
print(f"\nResults:")
print(f"  Successful traces: {stats['n_successful_traces']}/{len(start_points_re)}")
print(f"  Mean length: {stats['mean_length_re']:.2f} RE")
print(f"  Stop reasons: {stats['stop_reasons']}")

# Examine individual trace
if traces[0] is not None:
    trace = traces[0]
    combined = trace['combined']
    print(f"\nFirst trace details:")
    print(f"  Points: {combined['n_points']}")
    print(f"  Length: {combined['length']:.2f} RE")
    print(f"  Stop reasons: {combined['stop_reasons']}")
Example 2: Multi-Field Comparison
python
Copy code
import kamodo_ccmc.tools.vectorfieldtracer as vft

# Create tracers for different field types
tracers = {
    'magnetic': vft.create_magnetic_tracer(ko),
    'velocity': vft.create_velocity_tracer(ko),
    'current': vft.create_current_tracer(ko)
}

# Common starting points
start_points_re = [
    (3.0, 0.0, 0.0),    # Dayside
    (0.0, 3.0, 0.0),    # Dawn
    (-8.0, 0.0, 0.0),   # Nightside
    (0.0, -3.0, 0.0)    # Dusk
]

time_hours = 15.0
all_traces = {}

# Trace each field type
for field_name, tracer in tracers.items():
    print(f"\nTracing {field_name} field lines...")
    
    try:
        traces = tracer.trace_multiple_lines(
            start_points_re, 
            time_hours,
            max_steps=3000,
            direction='both'
        )
        
        stats = tracer.get_field_statistics(traces)
        all_traces[field_name] = traces
        
        print(f"  Success rate: {stats['n_successful_traces']}/{len(start_points_re)}")
        print(f"  Mean length: {stats['mean_length_re']:.2f} RE")
        print(f"  Field strength range: {stats['min_field_strength']:.2e} to {stats['max_field_strength']:.2e}")
        
    except Exception as e:
        print(f"  Error tracing {field_name}: {e}")

# Compare field line lengths
print(f"\nField line length comparison:")
for field_name, traces in all_traces.items():
    if traces:
        lengths = []
        for trace in traces:
            if trace and 'combined' in trace:
                lengths.append(trace['combined']['length'])
        
        if lengths:
            print(f"  {field_name.capitalize()}: {np.mean(lengths):.2f} ± {np.std(lengths):.2f} RE")
Example 3: Detailed Field Evaluation
python
Copy code
import kamodo_ccmc.tools.vectorfieldtracer as vft
import numpy as np

# Create tracer
tracer = vft.create_magnetic_tracer(ko)

# Test field evaluation at specific points
test_points_re = [
    (2.0, 0.0, 0.0),    # Close to Earth, dayside
    (5.0, 0.0, 0.0),    # Mid-distance, dayside
    (10.0, 0.0, 0.0),   # Far distance, dayside
    (-15.0, 0.0, 0.0),  # Nightside tail
    (0.0, 8.0, 2.0),    # Dawn sector, off-equatorial
]

time_hours = 12.0

print("Field evaluation at test points:")
print("Point (RE)          | Field Strength (nT) | Components (nT)")
print("-" * 65)

for point in test_points_re:
    try:
        # Get field components
        bx, by, bz = tracer.get_vector_field(point[0], point[1], point[2], time_hours)
        
        # Handle potential array returns
        if hasattr(bx, '__len__') and len(bx) == 1:
            bx, by, bz = float(bx[0]), float(by[0]), float(bz[0])
        
        b_mag = np.sqrt(bx**2 + by**2 + bz**2)
        
        print(f"({point[0]:5.1f}, {point[1]:5.1f}, {point[2]:5.1f}) | "
              f"{b_mag:15.2e} | ({bx:8.2e}, {by:8.2e}, {bz:8.2e})")
        
    except Exception as e:
        print(f"({point[0]:5.1f}, {point[1]:5.1f}, {point[2]:5.1f}) | "
              f"{'ERROR':15s} | {str(e)[:30]}")

# Trace from the most interesting point
best_point = test_points_re[2]  # 10 RE dayside
print(f"\nTracing field line from {best_point}:")

trace_result = tracer.trace_vector_line(
    best_point, time_hours,
    max_steps=5000,
    direction='both',
    verbose=True
)

if 'combined' in trace_result:
    combined = trace_result['combined']
    coords = combined['coordinates']
    
    print(f"  Total points: {len(coords)}")
    print(f"  Total length: {combined['length']:.2f} RE")
    print(f"  Start point: ({coords[0][0]:.2f}, {coords[0][1]:.2f}, {coords[0][2]:.2f}) RE")
    print(f"  End point: ({coords[-1][0]:.2f}, {coords[-1][1]:.2f}, {coords[-1][2]:.2f}) RE")
    print(f"  Stop reasons: {combined['stop_reasons']}")
Example 4: Custom Tracer Configuration
python
Copy code
import kamodo_ccmc.tools.vectorfieldtracer as vft

# Create custom tracer with specific settings
custom_tracer = vft.KamodoVectorFieldTracer(
    kamodo_object=ko,
    vector_components=('B_x', 'B_y', 'B_z'),
    coord_system='GSM',
    field_type='magnetic',
    auto_bounds=False,  # Don't use automatic bounds
    verbose=True        # Enable detailed output
)

# Manually set custom parameters
custom_tracer.max_distance = 25.0          # Max distance: 25 RE
custom_tracer.default_step_size = 0.05     # Smaller step size: 0.05 RE
custom_tracer.weak_field_threshold = 1e-15 # More sensitive weak field detection

# Set custom coordinate bounds (x_min, x_max, y_min, y_max, z_min, z_max)
custom_tracer.default_bounds = (-30, 15, -20, 20, -20, 20)

print("Custom tracer configuration:")
print(f"  Max distance: {custom_tracer.max_distance} RE")
print(f"  Step size: {custom_tracer.default_step_size} RE")
print(f"  Coordinate bounds: {custom_tracer.default_bounds}")

# Use custom tracer
start_point = (8.0, 0.0, 0.0)
trace_result = custom_tracer.trace_vector_line(
    start_point, 
    time_hours=12.0,
    max_steps=10000,
    adaptive_step=True
)

if 'combined' in trace_result:
    print(f"\nCustom trace results:")
    print(f"  Length: {trace_result['combined']['length']:.3f} RE")
    print(f"  Points: {trace_result['combined']['n_points']}")
    print(f"  Stop reasons: {trace_result['combined']['stop_reasons']}")
Example 5: Batch Processing with Error Handling
python
Copy code
import kamodo_ccmc.tools.vectorfieldtracer as vft
import numpy as np

def process_time_series(ko, time_points, start_points_re):
    """Process multiple time points with error handling."""
    
    # Create tracer
    tracer = vft.create_magnetic_tracer(ko)
    
    results = {}
    
    for time_hours in time_points:
        print(f"\nProcessing time {time_hours:.1f} hours...")
        
        try:
            # Trace field lines at this time
            traces = tracer.trace_multiple_lines(
                start_points_re, 
                time_hours,
                max_steps=3000,
                direction='both'
            )
            
            # Get statistics
            stats = tracer.get_field_statistics(traces)
            
            # Store results
            results[time_hours] = {
                'traces': traces,
                'stats': stats,
                'success': True
            }
            
            print(f"  Success: {stats['n_successful_traces']}/{len(start_points_re)} traces")
            print(f"  Mean length: {stats['mean_length_re']:.2f} RE")
            
        except Exception as e:
            print(f"  Error at time {time_hours}: {e}")
            results[time_hours] = {
                'traces': None,
                'stats': None,
                'success': False,
                'error': str(e)
            }
    
    return results

# Example usage
time_points = [6.0, 12.0, 18.0, 24.0]  # Four times of day
start_points_re = [
    (4.0, 0.0, 0.0),   # Dayside
    (-10.0, 0.0, 0.0), # Nightside  
    (0.0, 6.0, 0.0),   # Dawn
    (0.0, -6.0, 0.0)   # Dusk
]

# Process time series
time_series_results = process_time_series(ko, time_points, start_points_re)

# Analyze temporal evolution
print(f"\nTemporal Analysis:")
print("Time (h) | Success Rate | Mean Length (RE)")
print("-" * 40)

for time_hours in sorted(time_series_results.keys()):
    result = time_series_results[time_hours]
    
    if result['success']:
        stats = result['stats']
        success_rate = stats['n_successful_traces'] / len(start_points_re) * 100
        mean_length = stats['mean_length_re']
        print(f"{time_hours:7.1f} | {success_rate:11.1f}% | {mean_length:13.2f}")
    else:
        print(f"{time_hours:7.1f} | {'ERROR':11s} | {'N/A':13s}")
Example 6: Model Diagnostics and Troubleshooting
python
Copy code
import kamodo_ccmc.tools.vectorfieldtracer as vft

# Step 1: Inspect your Kamodo object
print("=== KAMODO OBJECT INSPECTION ===")
available_vars = vft.inspect_kamodo_object(ko)

# Step 2: Test basic integration
print("\n=== INTEGRATION TEST ===")
vft.test_kamodo_integration(ko, test_time=12.0)

# Step 3: Get model bounds
print("\n=== MODEL BOUNDS ===")
if 'B_x' in available_vars:
    bounds = vft.get_model_bounds(ko, ('B_x', 'B_y', 'B_z'))
    print(f"Model boundaries:")
    print(f"  X: {bounds['x_range'][0]:.1f} to {bounds['x_range'][1]:.1f} RE")
    print(f"  Y: {bounds['y_range'][0]:.1f} to {bounds['y_range'][1]:.1f} RE")
    print(f"  Z: {bounds['z_range'][0]:.1f} to {bounds['z_range'][1]:.1f} RE")
    print(f"  Max extent: {bounds['max_extent']:.1f} RE")

# Step 4: Run complete example
print("\n=== COMPLETE EXAMPLE ===")
try:
    results = vft.example_usage_with_kamodo(ko)
    
    if results:
        print("✅ Example completed successfully!")
        print(f"Created tracers: {list(results['tracers'].keys())}")
        print(f"Total traces: {sum(len(traces) for traces in results['traces'].values())}")
    else:
        print("❌ Example failed - check your Kamodo object")
        
except Exception as e:
    print(f"❌ Example error: {e}")
    print("Check that your Kamodo object contains the required vector field components")
Performance Tips
Use appropriate step sizes: Smaller steps = higher accuracy but slower execution
Set reasonable max_steps: Prevent infinite loops in complex field regions
Use adaptive stepping: Automatically adjusts step size based on field strength
Limit search domains: Use model bounds to avoid extrapolation errors
Batch processing: Trace multiple lines together for better efficiency
Troubleshooting
Common Issues
Issue: ValueError: Vector component 'B_x' not found in Kamodo object

Solution: Check available variables with vft.inspect_kamodo_object(ko)
Use the correct component names for your model
Issue: Field lines stop immediately with "Weak_field" reason

Solution: Adjust weak_field_threshold or check your starting points
Some models have different field strength scales
Issue: "Field_evaluation_error" during tracing

Solution: Ensure starting points are within model domain
Use vft.get_model_bounds() to check valid coordinate ranges
Issue: Traces hit "Model_boundary" quickly

Solution: Use auto_bounds=True or manually set appropriate bounds
Check that your model covers the regions you want to trace
Getting Help
Run diagnostics: Use vft.test_kamodo_integration(ko) to identify issues
Check examples: Run vft.example_usage_with_kamodo(ko) for a complete test
Inspect your model: Use vft.inspect_kamodo_object(ko) to understand structure
Enable verbose output: Set verbose=True to see detailed tracing information
Integration with Other Modules
The vector field tracer integrates seamlessly with other Kamodo tools:

python
Copy code
# Use with last closed field line finder
import kamodo_ccmc.tools.lastclosedmag as lcm
results = lcm.find_last_closed_field_lines(ko, time_hours=12.0)

# Use with null point finder  
import kamodo_ccmc.tools.nullfinder as nf
nulls, analysis = nf.find_magnetic_null_points(results, ko, time_hours=12.0)
This creates a complete workflow for magnetospheric analysis from basic field line tracing to advanced reconnection site detection.

