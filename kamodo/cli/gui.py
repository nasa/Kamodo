import hydra
from hydra.experimental import compose
from omegaconf import OmegaConf
import dash
import dash_core_components as dcc
import dash_html_components as html
from kamodo.cli.main import eval_config
from kamodo import get_defaults
import plotly.graph_objs as go
from dash.dependencies import Input, Output
import numpy as np
from os import path
import os
import numpy as np
from sympy.core.function import UndefinedFunction 
import pandas as pd
import sys 
from dash_katex import DashKatex

try:
    hydra.experimental.initialize(strict=False)
except:
    pass

from ast import literal_eval


def get_equation_divs(model, model_name, model_params):

    equations = []

    dcc.Checklist(
            options=[
                {'label': "NYC", 'value': 'NYC'},
                {'label': 'Montréal', 'value': 'MTL'},
                {'label': 'San Francisco', 'value': 'SF'}
            ],
            value=['MTL', 'SF']
        )
    options = []
    selected = []
    for var_symbol, func in model.items():
        if type(var_symbol) != UndefinedFunction:
            var_label = str(type(var_symbol))

            # make sure user-specified params get a check mark
            if str(type(var_symbol)) in model_params:
                selected.append(var_label)

            options.append({'label': var_label, 'value': var_label})
            equations.append(
                DashKatex(
                    id = "{}-{}-expression".format(model_name, str(var_symbol)),
                    expression = "{}".format(model.to_latex(
                        keys = [str(var_symbol)],
                        mode = "plain")),
                    )
                )

    equation_divs = html.Div([
            dcc.Checklist(
                options = options,
                value = selected,
                className = "six columns"
                ),
            html.Div(
                children = equations,
                className = "six columns"
                ),],
        className = 'row',
        )

    return equation_divs

def get_gui(cfg):

    # moved into assets
    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
    external_scripts = ['https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-MML-AM_CHTML']


    app = dash.Dash(__name__, 
        external_stylesheets = external_stylesheets,
        external_scripts = external_scripts,
        )


    app.layout = html.Div(children=[
        html.H1(children='Kamodo'),

        html.Div(children='''
            Kamodo: A low-coding interface for models and data
        '''),
    ])

    tabs_children = []
    graphs = {}

    for model_name, model_conf in cfg.models.items():
        print(model_name)

        try:
            model = hydra.utils.instantiate(model_conf)
        except Exception as m:
            print('could not instantiate {}'.format(model_name))
            print(m)

        if len(model) == 0:
            continue

        model_params = {}
        if model_conf.plot is None:
            for var_symbol, func in model.items():
                if type(var_symbol) == UndefinedFunction:
                    model_params[str(var_symbol)] = get_defaults(func)
        else:
            for var_name, var_args in model_conf.plot.items():
                if var_args is None:
                    try:
                        model_params[var_name] = get_defaults(model[var_name])
                    except Exception as m:
                        print(model.detail)
                        print(m)
                        sys.exit()
                else:
                    model_params[var_name] = eval_config(var_args)


        print(model_params)
     

        tab_children = []

        tab_children.append(get_equation_divs(model, model_name, model_params))

        
        for varname, params in model_params.items():
            # plot_args = {varname: eval_config(params)}
            plot_args = {varname: params}
            try:
                fig = model.figure(
                    varname,
                    indexing = 'ij',
                    return_type = True, 
                    **plot_args[varname])
                
                chart_type = fig.get('chart_type', None)
                figure = go.Figure(fig)
                if chart_type == 'line':
                    axis_type = type(list(plot_args[varname].values())[0])
                    if axis_type == pd.DatetimeIndex:
                        print('found datetime!')
                        figure.update_layout(
                            xaxis=go.layout.XAxis(
                                rangeslider=dict(
                                    visible=True
                                ),
                                type="date"
                            )
                        )
                    else:
                        print("asix type:", axis_type, plot_args)
                figure.layout['width'] = None
                figure.layout['height'] = None
                figure.update_layout(autosize = False,
                    # width = '100%', height='100%',
                    )

                graph_id = "graph-{}-{}".format(model_name, varname)
                
                print('storing graph')
                graphs[graph_id] = [
                    chart_type,
                    varname,
                    plot_args[varname],
                    model,
                    dcc.Graph(
                        id = graph_id,
                        figure = figure,
                        style = dict(height = 'inherit', width = '800px'))]

                print("len(graphs) {}".format(len(graphs)))

                tab_children.append(
                    html.Div(
                        html.Div(
                            graphs[graph_id][-1],
                            className = 'twelve columns'),
                        className = 'row'))

                # if chart_type == 'line':
                #     tab_children.extend(get_range_slider(
                #         graph_id, **plot_args[varname]))


            except Exception as m:
                print('could not plot {} with params:'.format(varname))
                print(plot_args)
                raise
                # print(m)


        tab = dcc.Tab(label = model_name, children = tab_children)
        
        tabs_children.append(tab)

    tabs = dcc.Tabs(children = tabs_children)
    app.layout.children.append(tabs)

    # generate_callbacks(app, graphs)
    
    return app

def generate_callbacks(app, graphs):
    print("generating callbacks for {} graphs".format(len(graphs)))
    for graph_id, (chart_type, varname, plot_args, model, graph) in graphs.items():
        if chart_type == 'line':
            generate_range_slider_callback(
                app,
                graph_id,
                varname,
                model,
                **plot_args)
        else:
            print("cannot generate callback for chart_type {}".format(chart_type))

def generate_range_slider_callback(app, graph_id, varname, model, **plot_kwargs):
    print("generating calback for {}: {}".format(graph_id, varname))
    output = Output(graph_id, "figure")
    input_dict = dict()
    for k,v in plot_kwargs.items():
        input_dict[k] = [
            Input("rangeslider-{}-{}".format(graph_id, k), "value"),
            Input("input-{}-{}".format(graph_id, k), "value")]

    inputs  = []
    for inputs_ in input_dict.values():
        inputs.extend(inputs_)

    @app.callback(output, inputs)
    def update_figure(*inputs):
        plot_args = dict()
        for i,k in enumerate(input_dict):
            plot_args[k] = np.linspace(*inputs[2*i], num = inputs[2*i+1])

        fig = model.figure(
                varname,
                indexing = 'ij',
                return_type = True, 
                **plot_args)
        return go.Figure(fig)


def get_range_slider(graph_id, **kwargs):
    divs = []
    for k,v in kwargs.items():
        try:
            print('type is: {}'.format(type(v[0])))
            marks = {i: "{0:0.2f}".format(i) for i in v[::int(len(v)/31)]}
        except:
            marks = {i: str(i) for i in v[::int(len(v)/31)]}

        # print(len(marks))
        # print(marks)

        row_children = [
            html.Div(
                className = "eight columns",
                children = [
                    dcc.RangeSlider(
                        id = "rangeslider-{}-{}".format(graph_id, k),
                        marks=marks,
                        min=min(v),
                        max=max(v),
                        step = (max(v) - min(v))/len(v),
                        value=[min(v), max(v)],
                        updatemode='drag',
                        )]),
            html.Div(
                className = "four columns",
                children = [
                    dcc.Input(
                        id = "input-{}-{}".format(graph_id, k),
                        placeholder='number of points in {}'.format(k),
                        type='number',
                        value = len(v),
                        min = 3,
                        max = 10*len(v))])
            ]

        row_div = html.Div(
            className = 'row',
            children = row_children)
        divs.append(row_div)

    return divs


def config_override(cfg):
    """Overrides with user-supplied configuration

    hourly will override its configuration using
    hourly.yaml if it is in the current working directory
    or users can set an override config:
        config_override=path/to/myconfig.yaml
    """
    if cfg.config_override is not None:
        override_path = hydra.utils.to_absolute_path(cfg.config_override)
        if path.exists(override_path):
            override_conf = OmegaConf.load(override_path)
            # merge overrides first input with second
            cfg = OmegaConf.merge(cfg, override_conf)
    return cfg

def main():

    cfg = compose('conf/kamodo.yaml')
    
    config_override = None
    if cfg.config_override is not None:
        override_path = "{}/{}".format(os.getcwd(),cfg.config_override)
        if path.exists(override_path):
            print("found {}".format(override_path))
            config_override = OmegaConf.load(override_path)
            print(config_override.pretty())
        else:
            if cfg.verbose > 0:
                print("could not get override: {}".format(override_path))
    

        if config_override is not None:
            cfg = OmegaConf.merge(cfg, config_override)
    

    if cfg.verbose > 0:
        print(cfg.pretty())

    app = get_gui(cfg)
    app.run_server(debug = True)

# entrypoint for package installer
def entry():
    main()


if __name__ == "__main__":
    main()
