"""
Kamodo Vector Field Tracer - Core field line tracing functionality

High-performance vector field tracer for Kamodo models.
Supports magnetic field lines, velocity streamlines, current streamlines, etc.
All coordinates are in Earth radii (RE).
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
from scipy.spatial.distance import cdist


def get_colormap(name):
    """
    Get colormap in a way that's compatible with different matplotlib versions.

    Parameters:
    -----------
    name : str
        Name of the colormap

    Returns:
    --------
    cmap : matplotlib colormap
        The requested colormap
    """
    try:
        # Try new interface first (matplotlib >= 3.7)
        import matplotlib
        if hasattr(matplotlib, 'colormaps'):
            return matplotlib.colormaps[name]
        else:
            # Fallback for matplotlib < 3.7
            return plt.cm.get_cmap(name)
    except (AttributeError, KeyError):
        # Additional fallback
        try:
            return plt.get_cmap(name)
        except:
            # Last resort - return a default colormap
            print(f"Warning: Could not find colormap '{name}', using 'viridis'")
            return plt.get_cmap('viridis')


def plot_vector_field_traces(traces, tracer, title_suffix="", backend='matplotlib',
                            show_earth=True, interactive=True):
    """
    Plot vector field traces with field-type specific styling.

    Parameters:
    -----------
    traces : list
        List of trace result dictionaries
    tracer : KamodoVectorFieldTracer
        The tracer object used to generate the traces
    title_suffix : str
        Additional text to add to plot titles
    backend : str
        'matplotlib' or 'plotly' for plotting backend
    show_earth : bool
        Whether to show Earth sphere/circle
    interactive : bool
        For plotly: whether to show interactive features

    Returns:
    --------
    fig : matplotlib.figure.Figure or plotly.graph_objects.Figure
        The generated figure
    """
    if backend.lower() == 'plotly':
        return _plot_vector_field_traces_plotly(
            traces, tracer, title_suffix, show_earth, interactive)
    else:
        return _plot_vector_field_traces_matplotlib(
            traces, tracer, title_suffix, show_earth)


def _plot_vector_field_traces_matplotlib(traces, tracer, title_suffix="", show_earth=True):
    """
    Create matplotlib version of vector field trace plots with aspect ratio 1.
    """
    fig = plt.figure(figsize=(16, 12))

    # Create subplots
    ax1 = fig.add_subplot(221, projection='3d')  # 3D view
    ax2 = fig.add_subplot(222)  # XY plane
    ax3 = fig.add_subplot(223)  # XZ plane
    ax4 = fig.add_subplot(224)  # Field strength along trace

    # Field-type specific colors and labels
    field_colors = {
        'magnetic': 'plasma',
        'velocity': 'viridis',
        'current': 'inferno',
        'electric': 'cividis'
    }

    field_units = {
        'magnetic': 'nT',
        'velocity': 'km/s',
        'current': 'A/m²',
        'electric': 'V/m'
    }

    colormap = field_colors.get(tracer.field_type, 'plasma')
    field_unit = field_units.get(tracer.field_type, 'units')

    # Use updated colormap access
    cmap = get_colormap(colormap)
    colors = cmap(np.linspace(0, 1, len(traces)))

    all_field_strengths = []
    all_coords = []  # For aspect ratio calculation

    for i, trace in enumerate(traces):
        if trace is None:
            continue

        # Use combined trace if available
        trace_data = trace.get('combined', trace.get('forward', trace.get('backward')))
        if trace_data is None or len(trace_data['coordinates']) < 2:
            continue

        coords = trace_data['coordinates']  # Already in RE
        field_mag = trace_data['field_magnitude']
        all_coords.extend(coords)

        # Coordinates are already in RE
        x_re, y_re, z_re = coords.T

        # 3D plot colored by field strength
        scatter = ax1.scatter(x_re, y_re, z_re, c=field_mag,
                              cmap=colormap, s=2, alpha=0.7)
        ax1.plot(x_re, y_re, z_re, color=colors[i], alpha=0.6, linewidth=1)

        # 2D projections
        ax2.plot(x_re, y_re, color=colors[i], alpha=0.8, linewidth=1.5)
        ax3.plot(x_re, z_re, color=colors[i], alpha=0.8, linewidth=1.5)

        # Field strength along trace
        s_coord = np.arange(len(field_mag)) * tracer.default_step_size  # Already in RE
        ax4.plot(s_coord, field_mag, color=colors[i], alpha=0.8, linewidth=1.5)

        all_field_strengths.extend(field_mag)

    # Add Earth
    if show_earth:
        theta = np.linspace(0, 2*np.pi, 100)
        for ax in [ax2, ax3]:
            ax.plot(np.cos(theta), np.sin(theta), 'b-', linewidth=2, alpha=0.7, label='Earth')
            ax.set_aspect('equal')  # Aspect ratio 1
            ax.grid(True, alpha=0.3)

        # 3D Earth
        u = np.linspace(0, 2 * np.pi, 30)
        v = np.linspace(0, np.pi, 30)
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

    # Labels and titles
    coord_sys = tracer.coord_system
    field_name = tracer.field_type.title()

    ax1.set_xlabel(f'X ({coord_sys}) [RE]')
    ax1.set_ylabel(f'Y ({coord_sys}) [RE]')
    ax1.set_zlabel(f'Z ({coord_sys}) [RE]')
    ax1.set_title(f'3D {field_name} Field Lines {title_suffix}')

    ax2.set_xlabel(f'X ({coord_sys}) [RE]')
    ax2.set_ylabel(f'Y ({coord_sys}) [RE]')
    ax2.set_title('Equatorial Plane (Z=0)')

    ax3.set_xlabel(f'X ({coord_sys}) [RE]')
    ax3.set_ylabel(f'Z ({coord_sys}) [RE]')
    ax3.set_title('Meridional Plane (Y=0)')

    ax4.set_xlabel('Distance along trace [RE]')
    ax4.set_ylabel(f'{field_name} Field Strength [{field_unit}]')
    ax4.set_title(f'{field_name} Field Strength vs Distance')
    ax4.grid(True, alpha=0.3)

    # Colorbar for 3D plot
    if all_field_strengths:
        cbar = plt.colorbar(scatter, ax=ax1, shrink=0.8)
        cbar.set_label(f'{field_name} Field Strength [{field_unit}]')

    plt.tight_layout()
    return fig


def _plot_vector_field_traces_plotly(traces, tracer, title_suffix="",
                                    show_earth=True, interactive=True):
    """
    Create Plotly version of vector field trace plots with aspect ratio 1.
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        import plotly.express as px
    except ImportError:
        raise ImportError("Plotly is required for plotly backend. Install with: pip install plotly")

    # Field-type specific colors and labels
    field_colors = {
        'magnetic': 'Plasma',
        'velocity': 'Viridis',
        'current': 'Inferno',
        'electric': 'Cividis'
    }

    field_units = {
        'magnetic': 'nT',
        'velocity': 'km/s',
        'current': 'A/m²',
        'electric': 'V/m'
    }

    colorscale = field_colors.get(tracer.field_type, 'Plasma')
    field_unit = field_units.get(tracer.field_type, 'units')
    coord_sys = tracer.coord_system
    field_name = tracer.field_type.title()

    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{"type": "scene", "rowspan": 2}, {"type": "xy"}],
               [None, {"type": "xy"}]],
        subplot_titles=[
            f'3D {field_name} Field Lines {title_suffix}',
            f'Equatorial Plane (Z=0)',
            f'Meridional Plane (Y=0)'
        ],
        horizontal_spacing=0.1,
        vertical_spacing=0.1
    )

    # Process traces
    all_field_strengths = []
    trace_colors = px.colors.sample_colorscale(colorscale, len(traces))

    for i, trace in enumerate(traces):
        if trace is None:
            continue

        # Use combined trace if available
        trace_data = trace.get('combined', trace.get('forward', trace.get('backward')))
        if trace_data is None or len(trace_data['coordinates']) < 2:
            continue

        coords = trace_data['coordinates']  # Already in RE
        field_mag = trace_data['field_magnitude']

        # Coordinates are already in RE
        x_re, y_re, z_re = coords.T

        # Create hover text
        hover_text = [
            f'Point {j}<br>'
            f'X: {x_re[j]:.2f} RE<br>'
            f'Y: {y_re[j]:.2f} RE<br>'
            f'Z: {z_re[j]:.2f} RE<br>'
            f'{field_name}: {field_mag[j]:.2e} {field_unit}'
            for j in range(len(x_re))
        ]

        # 3D trace
        fig.add_trace(
            go.Scatter3d(
                x=x_re, y=y_re, z=z_re,
                mode='lines+markers',
                marker=dict(
                    size=2,
                    color=field_mag,
                    colorscale=colorscale,
                    opacity=0.8,
                    colorbar=dict(
                        title=f'{field_name}<br>[{field_unit}]',
                        x=0.45, len=0.5
                    ) if i == 0 else None,
                    showscale=i == 0
                ),
                line=dict(color=trace_colors[i], width=3),
                name=f'Trace {i+1}',
                hovertext=hover_text,
                hoverinfo='text',
                showlegend=False
            ),
            row=1, col=1
        )

        # XY projection (equatorial plane)
        fig.add_trace(
            go.Scatter(
                x=x_re, y=y_re,
                mode='lines',
                line=dict(color=trace_colors[i], width=2),
                name=f'Trace {i+1}',
                hovertemplate=f'X: %{{x:.2f}} RE<br>Y: %{{y:.2f}} RE<extra></extra>',
                showlegend=False
            ),
            row=1, col=2
        )

        # XZ projection (meridional plane)
        fig.add_trace(
            go.Scatter(
                x=x_re, y=z_re,
                mode='lines',
                line=dict(color=trace_colors[i], width=2),
                name=f'Trace {i+1}',
                hovertemplate=f'X: %{{x:.2f}} RE<br>Z: %{{y:.2f}} RE<extra></extra>',
                showlegend=False
            ),
            row=2, col=2
        )

        all_field_strengths.extend(field_mag)

    # Add Earth if requested
    if show_earth:
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

        # 2D Earth circles
        theta = np.linspace(0, 2*np.pi, 100)
        earth_x = np.cos(theta)
        earth_y = np.sin(theta)

        # XY plane Earth
        fig.add_trace(
            go.Scatter(
                x=earth_x, y=earth_y,
                mode='lines',
                line=dict(color='blue', width=3),
                name='Earth',
                hoverinfo='skip',
                showlegend=False
            ),
            row=1, col=2
        )

        # XZ plane Earth
        fig.add_trace(
            go.Scatter(
                x=earth_x, y=earth_y,
                mode='lines',
                line=dict(color='blue', width=3),
                name='Earth',
                hoverinfo='skip',
                showlegend=False
            ),
            row=2, col=2
        )

    # Update layout with aspect ratio 1
    fig.update_layout(
        title=f'{field_name} Field Line Traces {title_suffix}',
        showlegend=False,
        height=800,
        width=1200,
        scene_aspectmode='data'  # This sets aspect ratio to 1 for 3D plot
    )

    # Update 3D scene
    fig.update_scenes(
        xaxis_title=f'X ({coord_sys}) [RE]',
        yaxis_title=f'Y ({coord_sys}) [RE]',
        zaxis_title=f'Z ({coord_sys}) [RE]',
        camera=dict(
            eye=dict(x=1.5, y=1.5, z=1.5)
        )
    )

    # Update 2D plots with aspect ratio 1
    fig.update_xaxes(title_text=f'X ({coord_sys}) [RE]', row=1, col=2)
    fig.update_yaxes(title_text=f'Y ({coord_sys}) [RE]', row=1, col=2,
                     scaleanchor="x", scaleratio=1)  # Aspect ratio 1

    fig.update_xaxes(title_text=f'X ({coord_sys}) [RE]', row=2, col=2)
    fig.update_yaxes(title_text=f'Z ({coord_sys}) [RE]', row=2, col=2,
                     scaleanchor="x2", scaleratio=1)  # Aspect ratio 1

    # Add grid
    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')

    if not interactive:
        fig.update_layout(
            dragmode=False,
            scrollZoom=False,
            doubleClick=False,
            showTips=False
        )

    return fig


def plot_field_strength_along_traces_plotly(traces, tracer, title_suffix=""):
    """
    Create a separate Plotly plot showing field strength along traces.

    Parameters:
    -----------
    traces : list
        List of trace result dictionaries
    tracer : KamodoVectorFieldTracer
        The tracer object used to generate the traces
    title_suffix : str
        Additional text to add to plot title

    Returns:
    --------
    fig : plotly.graph_objects.Figure
        Plotly figure showing field strength vs distance
    """
    try:
        import plotly.graph_objects as go
        import plotly.express as px
    except ImportError:
        raise ImportError("Plotly is required. Install with: pip install plotly")

    field_units = {
        'magnetic': 'nT',
        'velocity': 'km/s',
        'current': 'A/m²',
        'electric': 'V/m'
    }

    field_unit = field_units.get(tracer.field_type, 'units')
    field_name = tracer.field_type.title()

    fig = go.Figure()

    colors = px.colors.qualitative.Plotly

    for i, trace in enumerate(traces):
        if trace is None:
            continue

        # Use combined trace if available
        trace_data = trace.get('combined', trace.get('forward', trace.get('backward')))
        if trace_data is None or len(trace_data['coordinates']) < 2:
            continue

        field_mag = trace_data['field_magnitude']

        # Distance along trace
        s_coord = np.arange(len(field_mag)) * tracer.default_step_size

        fig.add_trace(
            go.Scatter(
                x=s_coord,
                y=field_mag,
                mode='lines',
                name=f'Trace {i+1}',
                line=dict(color=colors[i % len(colors)], width=2),
                hovertemplate=(f'Distance: %{{x:.2f}} RE<br>{field_name}: '
                              f'%{{y:.2e}} {field_unit}<extra></extra>')
            )
        )

    fig.update_layout(
        title=f'{field_name} Field Strength Along Traces {title_suffix}',
        xaxis_title='Distance Along Trace [RE]',
        yaxis_title=f'{field_name} Field Strength [{field_unit}]',
        hovermode='closest',
        showlegend=True,
        width=800,
        height=500
    )

    fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')

    return fig


def plot_traces_interactive(traces, tracer, title_suffix="", show_field_strength=True):
    """
    Create interactive Plotly plots for vector field traces.

    Parameters:
    -----------
    traces : list
        List of trace result dictionaries
    tracer : KamodoVectorFieldTracer
        The tracer object used to generate the traces
    title_suffix : str
        Additional text to add to plot titles
    show_field_strength : bool
        Whether to also create field strength plot

    Returns:
    --------
    figs : dict
        Dictionary containing 'main' figure and optionally 'field_strength' figure
    """
    figs = {}

    # Main trace plot
    figs['main'] = plot_vector_field_traces(
        traces, tracer, title_suffix=title_suffix,
        backend='plotly', interactive=True
    )

    # Field strength plot
    if show_field_strength:
        figs['field_strength'] = plot_field_strength_along_traces_plotly(
            traces, tracer, title_suffix=title_suffix
        )

    return figs

def get_model_bounds(ko, vector_components, verbose=True):
    """
    Get the coordinate bounds from Kamodo model for the given vector components.

    Parameters:
    -----------
    ko : Kamodo object
        Kamodo object containing vector field data
    vector_components : tuple of str
        Vector component names (e.g., ('B_x', 'B_y', 'B_z'))
    verbose : bool
        Whether to print bounds information

    Returns:
    --------
    bounds : dict
        Dictionary containing coordinate bounds
    """
    try:
        import kamodo_ccmc.flythrough.model_wrapper as MW

        # Use the first vector component to get coordinate ranges
        var = vector_components[0]

        if verbose:
            print(f"Getting coordinate bounds from Kamodo model using '{var}'...")

        cr = MW.Coord_Range(ko, [var], return_dict=True, print_output=False)

        xmin, xmax = cr[var]['X'][0], cr[var]['X'][1]
        ymin, ymax = cr[var]['Y'][0], cr[var]['Y'][1]
        zmin, zmax = cr[var]['Z'][0], cr[var]['Z'][1]

        bounds = {
            'x_range': (xmin, xmax),
            'y_range': (ymin, ymax),
            'z_range': (zmin, zmax),
            'x_center': (xmin + xmax) / 2,
            'y_center': (ymin + ymax) / 2,
            'z_center': (zmin + zmax) / 2,
            'max_extent': max(xmax - xmin, ymax - ymin, zmax - zmin) / 2
        }

        if verbose:
            print(f"Model coordinate bounds (in RE):")
            print(f"  X: [{xmin:.1f}, {xmax:.1f}] RE")
            print(f"  Y: [{ymin:.1f}, {ymax:.1f}] RE")
            print(f"  Z: [{zmin:.1f}, {zmax:.1f}] RE")
            print(f"  Max extent: {bounds['max_extent']:.1f} RE")

        return bounds

    except ImportError:
        print("Warning: kamodo_ccmc.flythrough.model_wrapper not available")
        print("Using default bounds")
        return {
            'x_range': (-60, 20),
            'y_range': (-40, 40),
            'z_range': (-40, 40),
            'x_center': -20,
            'y_center': 0,
            'z_center': 0,
            'max_extent': 40
        }
    except Exception as e:
        print(f"Error getting model bounds: {e}")
        print("Using default bounds")
        return {
            'x_range': (-60, 20),
            'y_range': (-40, 40),
            'z_range': (-40, 40),
            'x_center': -20,
            'y_center': 0,
            'z_center': 0,
            'max_extent': 40
        }


class KamodoVectorFieldTracer:
    """
    High-performance vector field tracer for Kamodo models.
    Supports magnetic field lines, velocity streamlines, current streamlines, etc.
    All coordinates are in Earth radii (RE).
    """

    def __init__(self, kamodo_object, vector_components, coord_system='GSM',
                 use_numba=True, earth_radius_km=6371.2, field_type='magnetic',
                 auto_bounds=True, verbose=False):
        """
        Initialize the vector field tracer for Kamodo models.

        Parameters:
        -----------
        kamodo_object : Kamodo
            Kamodo object (ko) containing vector field functions
        vector_components : tuple of str
            Names of vector components in Kamodo (e.g., ('B_x', 'B_y', 'B_z'))
        coord_system : str
            Coordinate system ('GSM', 'GSE', 'SM', etc.)
        use_numba : bool
            Whether to use Numba JIT compilation
        earth_radius_km : float
            Earth radius in km (default: 6371.2 km) - for reference only
        field_type : str
            Type of vector field ('magnetic', 'velocity', 'current', 'electric', 'other')
        auto_bounds : bool
            Whether to automatically determine bounds from model
        verbose : bool
            Whether to print detailed information
        """
        self.ko = kamodo_object
        self.vector_components = vector_components
        self.coord_system = coord_system
        self.use_numba = use_numba
        self.earth_radius_km = earth_radius_km  # For reference/conversion
        self.field_type = field_type.lower()
        self.verbose = verbose

        # Usage counter
        self.use_counter = 0

        # Validate vector components exist in Kamodo object
        available_vars = list(str(item) for item in kamodo_object.keys())
        for component in vector_components:
            if component not in available_vars:
                raise ValueError(f"Vector component '{component}' not found in Kamodo object")

        # Field-specific configuration (all in RE units)
        self._configure_field_parameters()

        # Get model bounds automatically if requested
        if auto_bounds:
            self.model_bounds = get_model_bounds(
                kamodo_object, vector_components, verbose=verbose)
            self._update_bounds_based_on_model()
        else:
            # Default magnetosphere-specific bounds (in RE)
            self.default_bounds = (-60, 20,  # X: -60 to 20 RE
                                   -40, 40,  # Y: -40 to 40 RE
                                   -40, 40)  # Z: -40 to 40 RE

    def _configure_field_parameters(self):
        """Configure parameters based on field type (all in RE units)."""
        if self.field_type == 'magnetic':
            self.min_altitude = 1.02  # 1.02 RE (about 120 km altitude)
            self.max_distance = 100   # 100 RE from Earth
            self.weak_field_threshold = 1e-12  # nT
            self.integration_method = 'field_line'  # Follow field lines
            self.default_step_size = 0.1  # 0.1 RE

        elif self.field_type == 'velocity':
            self.min_altitude = 1.02  # 1.02 RE
            self.max_distance = 50    # 50 RE (smaller for velocity)
            self.weak_field_threshold = 1e-3  # km/s
            self.integration_method = 'streamline'  # Follow streamlines
            self.default_step_size = 0.2  # 0.2 RE

        elif self.field_type == 'current':
            self.min_altitude = 1.02  # 1.02 RE
            self.max_distance = 30    # 30 RE (current systems are more localized)
            self.weak_field_threshold = 1e-9  # A/m²
            self.integration_method = 'streamline'
            self.default_step_size = 0.15  # 0.15 RE

        elif self.field_type == 'electric':
            self.min_altitude = 1.02  # 1.02 RE
            self.max_distance = 20    # 20 RE
            self.weak_field_threshold = 1e-6  # V/m
            self.integration_method = 'field_line'
            self.default_step_size = 0.1  # 0.1 RE

        else:  # 'other' or custom field
            self.min_altitude = 1.02  # 1.02 RE
            self.max_distance = 50    # 50 RE
            self.weak_field_threshold = 1e-12
            self.integration_method = 'streamline'
            self.default_step_size = 0.1  # 0.1 RE

    def _update_bounds_based_on_model(self):
        """Update tracing parameters based on model bounds."""
        if hasattr(self, 'model_bounds'):
            bounds = self.model_bounds

            # Update default bounds
            self.default_bounds = (bounds['x_range'][0], bounds['x_range'][1],
                                   bounds['y_range'][0], bounds['y_range'][1],
                                   bounds['z_range'][0], bounds['z_range'][1])

            # Update max distance (stay within model bounds)
            max_extent = bounds['max_extent']
            self.max_distance = min(self.max_distance, max_extent * 0.95)

            if self.verbose:
                print(f"Updated stopping criteria based on model bounds:")
                print(f"  Max distance: {self.max_distance:.1f} RE")

    def get_vector_field(self, x_re, y_re, z_re, time_hours):
        """
        Get vector field components at given position(s) and time using Kamodo format.

        Parameters:
        -----------
        x_re, y_re, z_re : float or array
            Position coordinates in Earth radii (RE)
        time_hours : float
            Time in hours since midnight of first day of data

        Returns:
        --------
        vx, vy, vz : float or array
            Vector field components
        """
        try:
            # Kamodo calling format: ko['component']([time, x, y, z])
            # Handle both single point and array inputs
            if np.isscalar(x_re) and np.isscalar(y_re) and np.isscalar(z_re):
                # Single point
                point = [time_hours, x_re, y_re, z_re]
                vx = self.ko[self.vector_components[0]](point)
                vy = self.ko[self.vector_components[1]](point)
                vz = self.ko[self.vector_components[2]](point)
            else:
                # Multiple points - ensure arrays
                x_re = np.atleast_1d(x_re)
                y_re = np.atleast_1d(y_re)
                z_re = np.atleast_1d(z_re)

                # Create points array with time prepended
                n_points = len(x_re)
                points = []
                for i in range(n_points):
                    points.append([time_hours, x_re[i], y_re[i], z_re[i]])

                # Get field components for all points
                vx = self.ko[self.vector_components[0]](*points)
                vy = self.ko[self.vector_components[1]](*points)
                vz = self.ko[self.vector_components[2]](*points)

        except Exception as e:
            raise ValueError(f"Could not evaluate vector field from Kamodo: {e}")

        return vx, vy, vz

    def get_vector_field_batch(self, positions_re, time_hours):
        """
        Get vector field components for multiple positions efficiently.

        Parameters:
        -----------
        positions_re : array_like
            Array of shape (n_points, 3) containing (x, y, z) coordinates in RE
        time_hours : float
            Time in hours since midnight of first day of data

        Returns:
        --------
        field_components : ndarray
            Array of shape (n_points, 3) containing (vx, vy, vz) components
        """
        positions_re = np.array(positions_re)
        n_points = len(positions_re)

        try:
            # Create points list for Kamodo with time prepended
            points = []
            for i in range(n_points):
                points.append([time_hours, positions_re[i, 0], positions_re[i, 1], positions_re[i, 2]])

            # Get all field components at once
            vx = self.ko[self.vector_components[0]](*points)
            vy = self.ko[self.vector_components[1]](*points)
            vz = self.ko[self.vector_components[2]](*points)

            # Ensure arrays
            vx = np.atleast_1d(vx)
            vy = np.atleast_1d(vy)
            vz = np.atleast_1d(vz)

            # Stack into (n_points, 3) array
            field_components = np.column_stack([vx, vy, vz])

        except Exception as e:
            raise ValueError(f"Could not evaluate batch vector field from Kamodo: {e}")

        return field_components

    def _check_stopping_criteria(self, x_re, y_re, z_re, field_magnitude):
        """
        Check if tracing should stop based on field type and position.

        Parameters:
        -----------
        x_re, y_re, z_re : float
            Current position in Earth radii (RE)
        field_magnitude : float
            Magnitude of vector field at current position

        Returns:
        --------
        should_stop : bool
            True if tracing should stop
        stop_reason : str
            Reason for stopping
        """
        r_re = np.sqrt(x_re**2 + y_re**2 + z_re**2)

        # Check if hit Earth's surface (with small buffer)
        if r_re < self.min_altitude:
            return True, "Earth_surface"

        # Check if too far from Earth
        if r_re > self.max_distance:
            return True, "Max_distance"

        # Check for weak field
        if field_magnitude < self.weak_field_threshold:
            return True, "Weak_field"

        # Field-specific stopping criteria
        if self.field_type == 'magnetic':
            # Magnetopause approximation
            if x_re > 0 and r_re > 15:  # Dayside (15 RE)
                return True, "Magnetopause_day"
            if x_re < -50:  # Nightside tail (50 RE)
                return True, "Magnetotail"

        elif self.field_type == 'velocity':
            # Velocity-specific boundaries (e.g., solar wind boundary)
            if r_re > 25:  # 25 RE
                return True, "Flow_boundary"

        elif self.field_type == 'current':
            # Current sheet boundaries
            if abs(z_re) < 0.5 and abs(x_re) > 10:  # Current sheet region
                return True, "Current_sheet"

        # Additional bounds checking based on model
        if hasattr(self, 'model_bounds'):
            bounds = self.model_bounds

            # Check if outside model bounds (with small margin)
            margin = 0.1  # 0.1 RE margin
            if (x_re < bounds['x_range'][0] + margin or
                x_re > bounds['x_range'][1] - margin or
                y_re < bounds['y_range'][0] + margin or
                y_re > bounds['y_range'][1] - margin or
                z_re < bounds['z_range'][0] + margin or
                z_re > bounds['z_range'][1] - margin):
                return True, "Model_boundary"

        return False, ""

    def trace_vector_line(self, start_point_re, time_hours, max_steps=10000,
                          step_size_re=None, adaptive_step=True, direction='both',
                          integration_method=None):
        """
        Trace a vector field line/streamline from starting point.

        Parameters:
        -----------
        start_point_re : tuple or array
            (x, y, z) starting coordinates in Earth radii (RE)
        time_hours : float
            Time in hours since midnight of first day of data
        max_steps : int
            Maximum number of integration steps
        step_size_re : float, optional
            Step size in RE (uses default if None)
        adaptive_step : bool
            Whether to use adaptive step sizing
        direction : str
            'forward', 'backward', or 'both'
        integration_method : str, optional
            Override default integration method

        Returns:
        --------
        trace_result : dict
            Dictionary containing trace results
        """
        self.use_counter += 1

        x0, y0, z0 = start_point_re

        if step_size_re is None:
            step_size_re = self.default_step_size

        if integration_method is None:
            integration_method = self.integration_method

        results = {}

        if direction in ['forward', 'both']:
            # Trace forward (along vector field)
            forward_trace = self._trace_single_direction(
                x0, y0, z0, time_hours, step_size_re, max_steps,
                adaptive_step, forward=True, method=integration_method)
            results['forward'] = forward_trace

        if direction in ['backward', 'both']:
            # Trace backward (against vector field)
            backward_trace = self._trace_single_direction(
                x0, y0, z0, time_hours, step_size_re, max_steps,
                adaptive_step, forward=False, method=integration_method)
            results['backward'] = backward_trace

        # Combine forward and backward if both traced
        if direction == 'both' and 'forward' in results and 'backward' in results:
            self._combine_traces(results)

        return results

    def _trace_single_direction(self, x0, y0, z0, time_hours, step_size_re,
                                max_steps, adaptive_step, forward=True, method='streamline'):
        """
        Trace vector field in single direction.

        Parameters:
        -----------
        x0, y0, z0 : float
            Starting coordinates in RE
        time_hours : float
            Time in hours
        step_size_re : float
            Step size in RE
        max_steps : int
            Maximum integration steps
        adaptive_step : bool
            Whether to use adaptive step sizing
        forward : bool
            Whether to trace in forward direction
        method : str
            Integration method ('field_line' or 'streamline')

        Returns:
        --------
        result : dict
            Dictionary with trace results
        """
        trajectory = np.zeros((max_steps + 1, 3))
        field_magnitude = np.zeros(max_steps + 1)
        field_components = np.zeros((max_steps + 1, 3))

        x, y, z = x0, y0, z0
        trajectory[0] = [x, y, z]

        # Initial field evaluation
        try:
            vx, vy, vz = self.get_vector_field(x, y, z, time_hours)
            # Ensure scalar values
            if np.isscalar(vx):
                v_mag = np.sqrt(vx**2 + vy**2 + vz**2)
                field_magnitude[0] = v_mag
                field_components[0] = [vx, vy, vz]
            else:
                # Handle case where Kamodo returns arrays even for single points
                vx, vy, vz = float(vx[0]), float(vy[0]), float(vz[0])
                v_mag = np.sqrt(vx**2 + vy**2 + vz**2)
                field_magnitude[0] = v_mag
                field_components[0] = [vx, vy, vz]
        except Exception as e:
            return self._create_trace_result(trajectory[:1], field_magnitude[:1],
                                             field_components[:1],
                                             f'Field_evaluation_error: {e}')

        current_step_size = step_size_re
        direction_factor = 1.0 if forward else -1.0
        total_length = 0.0

        for i in range(max_steps):
            # Check stopping criteria
            should_stop, stop_reason = self._check_stopping_criteria(x, y, z, v_mag)
            if should_stop:
                return self._create_trace_result(
                    trajectory[:i+1], field_magnitude[:i+1],
                    field_components[:i+1], stop_reason, total_length)

            # Get vector field
            try:
                vx, vy, vz = self.get_vector_field(x, y, z, time_hours)
                # Ensure scalar values
                if not np.isscalar(vx):
                    vx, vy, vz = float(vx[0]), float(vy[0]), float(vz[0])

                v_mag = np.sqrt(vx**2 + vy**2 + vz**2)
                field_magnitude[i+1] = v_mag
                field_components[i+1] = [vx, vy, vz]

            except Exception as e:
                return self._create_trace_result(
                    trajectory[:i+1], field_magnitude[:i+1],
                    field_components[:i+1],
                    f'Field_evaluation_error: {e}', total_length)

            if v_mag < self.weak_field_threshold:
                return self._create_trace_result(
                    trajectory[:i+1], field_magnitude[:i+1],
                    field_components[:i+1], 'Weak_field', total_length)

            # Adaptive step sizing
            if adaptive_step:
                current_step_size = self._compute_adaptive_step(
                    v_mag, step_size_re, method)

            # Integration step
            if method == 'field_line':
                # For field lines (e.g., magnetic, electric), follow field direction
                dx = direction_factor * current_step_size * vx / v_mag
                dy = direction_factor * current_step_size * vy / v_mag
                dz = direction_factor * current_step_size * vz / v_mag
            else:  # streamline
                # For streamlines (velocity, current), integrate as streamline
                dt = current_step_size / v_mag  # Time step
                dx = direction_factor * dt * vx
                dy = direction_factor * dt * vy
                dz = direction_factor * dt * vz

            # Update position
            x += dx
            y += dy
            z += dz
            trajectory[i+1] = [x, y, z]

            # Update total length
            step_length = np.sqrt(dx**2 + dy**2 + dz**2)
            total_length += step_length

        return self._create_trace_result(
            trajectory[:max_steps+1], field_magnitude[:max_steps+1],
            field_components[:max_steps+1], 'Max_steps', total_length)

    def _compute_adaptive_step(self, field_mag, base_step, method):
        """
        Compute adaptive step size based on field magnitude and type.

        Parameters:
        -----------
        field_mag : float
            Field magnitude
        base_step : float
            Base step size in RE
        method : str
            Integration method

        Returns:
        --------
        step_size : float
            Computed step size in RE
        """
        if self.field_type == 'magnetic':
            # Magnetic field adaptive stepping
            if field_mag > 1000:  # Strong field (nT)
                factor = 0.5
            elif field_mag < 10:  # Weak field
                factor = 2.0
            else:
                factor = 1.0
        elif self.field_type == 'velocity':
            # Velocity field adaptive stepping
            if field_mag > 500:  # High velocity (km/s)
                factor = 0.3
            elif field_mag < 50:  # Low velocity
                factor = 1.5
            else:
                factor = 1.0
        elif self.field_type == 'current':
            # Current field adaptive stepping
            if field_mag > 1e-6:  # Strong current
                factor = 0.4
            elif field_mag < 1e-8:  # Weak current
                factor = 2.0
            else:
                factor = 1.0
        else:
            factor = 1.0

        # Limit step size changes
        factor = np.clip(factor, 0.1, 3.0)
        return base_step * factor

    def _create_trace_result(self, trajectory, field_magnitude, field_components,
                             stop_reason, total_length=None):
        """
        Create standardized trace result dictionary.

        Parameters:
        -----------
        trajectory : ndarray
            Array of shape (n_points, 3) with coordinates
        field_magnitude : ndarray
            Array of shape (n_points,) with field magnitudes
        field_components : ndarray
            Array of shape (n_points, 3) with field components
        stop_reason : str
            Reason for stopping trace
        total_length : float, optional
            Total trajectory length in RE

        Returns:
        --------
        result : dict
            Trace result dictionary
        """
        if total_length is None:
            if len(trajectory) > 1:
                total_length = np.sum(np.linalg.norm(np.diff(trajectory, axis=0), axis=1))
            else:
                total_length = 0.0

        return {
            'coordinates': trajectory,  # In RE
            'field_magnitude': field_magnitude,
            'field_components': field_components,
            'stop_reason': stop_reason,
            'length': total_length,  # In RE
            'n_points': len(trajectory)
        }

    def _combine_traces(self, results):
        """
        Combine forward and backward traces.

        Parameters:
        -----------
        results : dict
            Dictionary with 'forward' and 'backward' trace results

        Updates the results dictionary with a 'combined' entry
        """
        if 'forward' not in results or 'backward' not in results:
            return

        # Get traces
        forward = results['forward']
        backward = results['backward']

        # Reverse backward trace and combine (remove duplicate starting point)
        back_coords = backward['coordinates'][::-1]
        forward_coords = forward['coordinates']

        if len(back_coords) > 0 and len(forward_coords) > 0:
            combined_coords = np.vstack([back_coords[:-1], forward_coords])

            back_field_mag = backward['field_magnitude'][::-1]
            forward_field_mag = forward['field_magnitude']
            combined_field_mag = np.concatenate([back_field_mag[:-1], forward_field_mag])

            back_field_comp = backward['field_components'][::-1]
            forward_field_comp = forward['field_components']
            combined_field_comp = np.vstack([back_field_comp[:-1], forward_field_comp])
        else:
            combined_coords = forward_coords if len(forward_coords) > 0 else back_coords
            combined_field_mag = (forward['field_magnitude'] if len(forward_coords) > 0
                                  else backward['field_magnitude'])
            combined_field_comp = (forward['field_components'] if len(forward_coords) > 0
                                  else backward['field_components'])

        # Calculate combined length
        combined_length = forward['length'] + backward['length']

        results['combined'] = {
            'coordinates': combined_coords,  # In RE
            'field_magnitude': combined_field_mag,
            'field_components': combined_field_comp,
            'length': combined_length,  # In RE
            'n_points': len(combined_coords),
            'stop_reasons': [backward['stop_reason'], forward['stop_reason']]
        }

    def trace_multiple_lines(self, start_points_re, time_hours, **kwargs):
        """
        Trace multiple vector field lines/streamlines.

        Parameters:
        -----------
        start_points_re : array_like
            Array of shape (n_lines, 3) containing starting points in RE
        time_hours : float
            Time in hours since midnight of first day of data
        **kwargs : dict
            Additional arguments passed to trace_vector_line

        Returns:
        --------
        traces : list
            List of trace result dictionaries
        """
        start_points_re = np.array(start_points_re)
        traces = []

        print(f"Tracing {len(start_points_re)} {self.field_type} field lines "
              f"at t={time_hours:.2f} hours...")

        for i, point in enumerate(start_points_re):
            try:
                trace_result = self.trace_vector_line(point, time_hours, **kwargs)
                traces.append(trace_result)

                # Progress update
                if (i + 1) % max(1, len(start_points_re) // 10) == 0:
                    print(f"  Completed {i + 1}/{len(start_points_re)} traces")

            except Exception as e:
                print(f"Error tracing from {point}: {e}")
                traces.append(None)

        return traces

    def get_field_statistics(self, trace_results):
        """
        Compute statistics for traced field lines.

        Parameters:
        -----------
        trace_results : list
            List of trace results

        Returns:
        --------
        stats : dict
            Statistics dictionary
        """
        lengths = []
        n_points = []
        max_field_strengths = []
        min_field_strengths = []
        stop_reasons = []

        for result in trace_results:
            if result is None:
                continue

            # Use combined trace if available
            trace_data = result.get('combined', result.get('forward', result.get('backward')))
            if trace_data is None:
                continue

            lengths.append(trace_data['length'])  # In RE
            n_points.append(trace_data['n_points'])

            field_mag = trace_data['field_magnitude']
            if len(field_mag) > 0:
                max_field_strengths.append(np.max(field_mag))
                min_field_strengths.append(np.min(field_mag[field_mag > 0]))

            if 'stop_reasons' in trace_data:
                stop_reasons.extend(trace_data['stop_reasons'])
            else:
                stop_reasons.append(trace_data['stop_reason'])

        # Count stop reasons
        unique_reasons, reason_counts = np.unique(stop_reasons, return_counts=True)
        reason_dict = dict(zip(unique_reasons, reason_counts))

        stats = {
            'n_successful_traces': len(lengths),
            'mean_length_re': np.mean(lengths) if lengths else 0,  # In RE
            'std_length_re': np.std(lengths) if lengths else 0,
            'mean_points': np.mean(n_points) if n_points else 0,
            'max_field_strength': np.max(max_field_strengths) if max_field_strengths else 0,
            'min_field_strength': np.min(min_field_strengths) if min_field_strengths else 0,
            'stop_reasons': reason_dict,
            'field_type': self.field_type,
            'vector_components': self.vector_components
        }

        return stats

    def reset_trace_usage(self):
        """
        Reset the number of times the tracer has been called to 0.

        Parameters:
        -----------

        Returns:
        --------
        """
        self.use_counter = 0
        print('Tracer usage reset to 0.')
 
        return

    def report_trace_usage(self, value=False):
        """
        Report the number of times the tracer has been called.

        Parameters:
        -----------
        value : bool
            Whether to return a value or to just print out usage.

        Returns:
        --------
        use_counter value or nothing, depending to value bool
        """
        if value:
            return self.use_counter

        print('Tracer usage: ',self.use_counter)

        return


# Convenience functions for creating tracers with correct component names
def create_magnetic_tracer(ko, coord_system='GSM', **kwargs):
    """Create magnetic field tracer with correct component names."""
    return KamodoVectorFieldTracer(
        kamodo_object=ko,
        vector_components=('B_x', 'B_y', 'B_z'),
        coord_system=coord_system,
        field_type='magnetic',
        **kwargs
    )


def create_velocity_tracer(ko, coord_system='GSM', **kwargs):
    """Create velocity field tracer with correct component names."""
    return KamodoVectorFieldTracer(
        kamodo_object=ko,
        vector_components=('v_x', 'v_y', 'v_z'),
        coord_system=coord_system,
        field_type='velocity',
        **kwargs
    )


def create_current_tracer(ko, coord_system='GSM', **kwargs):
    """Create current field tracer with correct component names."""
    return KamodoVectorFieldTracer(
        kamodo_object=ko,
        vector_components=('J_x', 'J_y', 'J_z'),
        coord_system=coord_system,
        field_type='current',
        **kwargs
    )


def create_electric_tracer(ko, coord_system='GSM', **kwargs):
    """Create electric field tracer with correct component names."""
    return KamodoVectorFieldTracer(
        kamodo_object=ko,
        vector_components=('E_x', 'E_y', 'E_z'),
        coord_system=coord_system,
        field_type='electric',
        **kwargs
    )


# Utility functions for diagnostics and setup
def inspect_kamodo_object(ko):
    """
    Inspect the structure and contents of a Kamodo object.
    
    Parameters:
    -----------
    ko : Kamodo object
        Your Kamodo object to inspect
    """
    print("Kamodo Object Inspection")
    print("=" * 25)
    
    # List all available variables
    vars_list = list(str(item) for item in ko.keys())
    print(f"Total variables: {len(vars_list)}")
    print(f"Variable names: {vars_list}")
    
    # Group by vector components
    vector_groups = {
        'Magnetic': ['B_x', 'B_y', 'B_z'],
        'Velocity': ['v_x', 'v_y', 'v_z', 'V_x', 'V_y', 'V_z'],
        'Current': ['J_x', 'J_y', 'J_z'],
        'Electric': ['E_x', 'E_y', 'E_z'],
    }
    
    print(f"\nVector field availability:")
    for field_name, components in vector_groups.items():
        available_comps = [comp for comp in components if comp in vars_list]
        if len(available_comps) == 3:
            print(f"✓ {field_name} field: {available_comps}")
        elif available_comps:
            print(f"⚠ {field_name} field (incomplete): {available_comps}")
        else:
            print(f"✗ {field_name} field: not available")
    
    # Check for other common variables
    other_vars = [var for var in vars_list if not any(var in group for group in vector_groups.values())]
    if other_vars:
        print(f"\nOther variables: {other_vars}")
    
    return vars_list


def test_kamodo_integration(ko, test_time=12.0):
    """
    Test Kamodo integration with detailed diagnostics.
    
    Parameters:
    -----------
    ko : Kamodo object
        Your Kamodo object
    test_time : float
        Test time in hours
    """
    print("Kamodo Integration Test")
    print("=" * 30)
    
    # Test basic Kamodo calls
    print("Testing basic Kamodo function calls...")
    
    # Test different calling formats (coordinates in RE)
    test_point = [test_time, 2.0, 0.0, 0.0]  # [time, x_RE, y_RE, z_RE]
    
    available_vars = list(str(item) for item in ko.keys())
    print(f"Available variables: {available_vars}")
    
    # Test magnetic field components if available
    if 'B_x' in available_vars:
        try:
            print(f"\nTesting B_x at {test_point} (coordinates in RE):")
            bx_val = ko['B_x'](test_point)
            print(f"B_x value: {bx_val}")
            print(f"B_x type: {type(bx_val)}")
            
            # Test multiple points - FIX: Use list of lists instead of *args
            test_points = [
                [test_time, 2.0, 0.0, 0.0],  # 2 RE from center
                [test_time, 3.0, 0.0, 0.0]   # 3 RE from center
            ]
            print(f"\nTesting B_x at multiple points (coordinates in RE):")
            # Try different calling methods for multiple points
            try:
                # Method 1: Pass as separate arguments
                bx_multi = ko['B_x'](*test_points)
                print(f"B_x multiple values (method 1): {bx_multi}")
            except Exception as e1:
                try:
                    # Method 2: Call individually
                    bx_multi = [ko['B_x'](point) for point in test_points]
                    print(f"B_x multiple values (method 2): {bx_multi}")
                except Exception as e2:
                    print(f"Error with multiple points - Method 1: {e1}")
                    print(f"Error with multiple points - Method 2: {e2}")
            
        except Exception as e:
            print(f"Error testing B_x: {e}")
    
    # Test the tracer creation and basic functionality
    try:
        if all(comp in available_vars for comp in ['B_x', 'B_y', 'B_z']):
            print(f"\nTesting magnetic field tracer creation:")
            mag_tracer = create_magnetic_tracer(ko)
            
            # Test single point field evaluation (coordinates in RE)
            test_coords = (2.0, 0.0, 0.0)  # 2 RE from center
            print(f"Testing field evaluation at {test_coords} RE:")
            
            bx, by, bz = mag_tracer.get_vector_field(*test_coords, test_time)
            
            # FIX: Handle numpy array formatting issue
            if hasattr(bx, '__len__') and len(bx) == 1:
                bx, by, bz = float(bx[0]), float(by[0]), float(bz[0])
            elif hasattr(bx, '__len__'):
                bx, by, bz = float(bx), float(by), float(bz)
            
            b_mag = np.sqrt(bx**2 + by**2 + bz**2)
            
            print(f"Magnetic field components: Bx={bx:.3e}, By={by:.3e}, Bz={bz:.3e}")
            print(f"Magnetic field magnitude: {b_mag:.3e} nT")
            
            # Test short trace
            print(f"\nTesting short field line trace:")
            trace_result = mag_tracer.trace_vector_line(
                test_coords, test_time, 
                max_steps=10, 
                direction='forward'
            )
            
            if 'forward' in trace_result:
                forward_trace = trace_result['forward']
                print(f"Trace points: {forward_trace['n_points']}")
                print(f"Trace length: {forward_trace['length']:.2f} RE")
                print(f"Stop reason: {forward_trace['stop_reason']}")
            
    except Exception as e:
        print(f"Error testing tracer: {e}")
        import traceback
        traceback.print_exc()


def example_usage_with_kamodo(ko):
    """
    Complete example of using the tracer with your Kamodo object 'ko'.
    All coordinates are in Earth radii (RE).
    
    Parameters:
    -----------
    ko : Kamodo object
        Your Kamodo object containing vector field data
    """
    print("Kamodo Vector Field Tracer - Complete Usage Example")
    print("All coordinates in Earth radii (RE)")
    print("=" * 55)
    
    try:
        # First, inspect the Kamodo object
        print("Step 1: Inspecting Kamodo object...")
        available_vars = inspect_kamodo_object(ko)
        
        # Create tracers based on available components
        print("\nStep 2: Creating tracers...")
        tracers = {}
        
        if all(comp in available_vars for comp in ['B_x', 'B_y', 'B_z']):
            tracers['magnetic'] = create_magnetic_tracer(ko)
            print("✓ Magnetic field tracer created")
        
        if all(comp in available_vars for comp in ['v_x', 'v_y', 'v_z']):
            tracers['velocity'] = create_velocity_tracer(ko)
            print("✓ Velocity field tracer created")
        
        if all(comp in available_vars for comp in ['J_x', 'J_y', 'J_z']):
            tracers['current'] = create_current_tracer(ko)
            print("✓ Current field tracer created")
        
        if all(comp in available_vars for comp in ['E_x', 'E_y', 'E_z']):
            tracers['electric'] = create_electric_tracer(ko)
            print("✓ Electric field tracer created")
        
        if not tracers:
            print("⚠ No recognized vector field components found")
            return None
        
        # Set up tracing parameters (all in RE)
        print("\nStep 3: Setting up tracing parameters...")
        
        # Define starting points in different regions (in RE)
        start_points_re = [
            (3.0, 0.0, 0.0),          # Dayside, equatorial
            (0.0, 3.0, 0.0),          # Dawn flank
            (0.0, -3.0, 0.0),         # Dusk flank  
            (-10.0, 0.0, 0.0),        # Nightside
            (2.0, 2.0, 1.0),          # Off-equatorial, northern
            (2.0, -2.0, -1.0),        # Off-equatorial, southern
        ]
        
        time_hours = 12.5  # Time since midnight of first day
        
        print(f"Starting points: {len(start_points_re)} locations (in RE)")
        print(f"Time: {time_hours} hours since midnight")
        
        # Test field evaluation at starting points
        print("\nStep 4: Testing field evaluation...")
        for field_type, tracer in tracers.items():
            print(f"\n{field_type.capitalize()} field evaluation:")
            try:
                for i, point in enumerate(start_points_re[:2]):  # Test first 2 points
                    vx, vy, vz = tracer.get_vector_field(point[0], point[1], point[2], time_hours)
                    v_mag = np.sqrt(vx**2 + vy**2 + vz**2)
                    print(f"  Point {i+1} ({point[0]:.1f}, {point[1]:.1f}, {point[2]:.1f}) RE: |{field_type[0].upper()}| = {v_mag:.3e}")
            except Exception as e:
                print(f"  Error: {e}")
        
        # Perform actual tracing
        print("\nStep 5: Tracing field lines...")
        all_traces = {}
        
        for field_type, tracer in tracers.items():
            print(f"\nTracing {field_type} field lines...")
            try:
                traces = tracer.trace_multiple_lines(
                    start_points_re, 
                    time_hours,
                    max_steps=5000,
                    step_size_re=0.1,      # Step size in RE
                    direction='both',
                    adaptive_step=True
                )
                all_traces[field_type] = traces
                
                # Get statistics
                stats = tracer.get_field_statistics(traces)
                print(f"  Successfully traced: {stats['n_successful_traces']}/{len(start_points_re)}")
                print(f"  Mean length: {stats['mean_length_re']:.2f} RE")
                print(f"  Stop reasons: {stats['stop_reasons']}")
                
            except Exception as e:
                print(f"  Error tracing {field_type}: {e}")
        
        print(f"\nStep 6: Results summary")
        print(f"Created {len(tracers)} tracers")
        print(f"Traced {sum(len(traces) for traces in all_traces.values())} field lines")
        
        return {
            'tracers': tracers,
            'traces': all_traces,
            'start_points_re': start_points_re,
            'time_hours': time_hours
        }
        
    except Exception as e:
        print(f"Error in example: {e}")
        return None


if __name__ == "__main__":
    print("=" * 60)
    print("KAMODO VECTOR FIELD TRACER")
    print("High-Performance 3D Field Line Tracing for Magnetospheric Models")
    print("All coordinates in Earth radii (RE)")
    print("=" * 60)
    
    print("\n🚀 QUICK START GUIDE:")
    print("1. Load your Kamodo object: ko = your_kamodo_model")
    print("2. Create tracer: tracer = create_magnetic_tracer(ko)")  
    print("3. Trace lines: traces = tracer.trace_multiple_lines(start_points_re, time_hours)")
    
    print("\n🔧 DIAGNOSTIC FUNCTIONS:")
    print("• inspect_kamodo_object(ko) - Check available variables")
    print("• test_kamodo_integration(ko) - Test basic functionality") 
    print("• example_usage_with_kamodo(ko) - Complete example")
    
    print("\n📊 SUPPORTED VECTOR FIELDS:")
    print("• Magnetic Fields: B_x, B_y, B_z (field lines)")
    print("• Velocity Fields: v_x, v_y, v_z (streamlines)")
    print("• Current Fields: J_x, J_y, J_z (current streamlines)")
    print("• Electric Fields: E_x, E_y, E_z (field lines)")
    
    print("\n⚙️ KEY FEATURES:")
    print("• Adaptive step sizing for accuracy")
    print("• Magnetosphere-aware stopping criteria")
    print("• Forward/backward/bidirectional tracing")
    print("• High-performance numerical integration")
    print("• All coordinates in Earth radii (RE)")

