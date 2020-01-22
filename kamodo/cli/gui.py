import hydra
from hydra.experimental import compose
from omegaconf import OmegaConf
import dash
import dash_core_components as dcc
import dash_html_components as html
from kamodo.cli.main import eval_config
import plotly.graph_objs as go
from dash.dependencies import Input, Output
import numpy as np
from os import path
import os
import numpy as np

hydra.experimental.initialize(strict=False)


def get_gui(cfg):
    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

    external_scripts = ['https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-MML-AM_CHTML']


    app = dash.Dash(__name__, 
        external_stylesheets = external_stylesheets,
        external_scripts = external_scripts)


    app.layout = html.Div(children=[
        html.H1(children='Kamodo'),

        html.Div(children='''
            Kamodo: A low-coding interface for models and data?
        '''),
    ])

    tabs_children = []

    for model_name, model_conf in cfg.models.items():
        print(model_name)
        

        try:
            model = hydra.utils.instantiate(model_conf)
        except Exception as m:
            print('could not instantiate {}'.format(model_name))
            print(m)

     
        tab_children = []

        for varname, params in model_conf.plot.items():
            plot_args = {varname: eval_config(params)}
            try:
                fig = model.figure(
                    varname,
                    indexing = 'ij',
                    return_type = True, 
                    **plot_args[varname])
                
                # fig = go.Figure(model.plot(**plot_args))
                chart_type = fig.get('chart_type', None)
                figure = go.Figure(fig)
                figure.layout['width'] = None
                figure.layout['height'] = None
                figure.update_layout(autosize = False,
                    # width = '100%', height='100%',
                    )

                # print(getattr(fig,'chart_type', 'no chart type'))
                tab_children.append(
                    html.Div(
                        html.Div(
                            dcc.Graph(
                                id = "graph-{}-{}".format(model_name, varname),
                                figure = figure,
                                style = dict(height = 'inherit', width = '800px'),
                                ),
                            className = 'twelve columns'),
                        className = 'row'))

                if chart_type == 'line':
                    tab_children.extend(get_range_slider(
                        model_name, varname, **plot_args[varname]))
                    # generate_range_slider_callback(
                    #     app,
                    #     model,
                    #     model_name,
                    #     varname,
                    #     **plot_args[varname])

            except Exception as m:
                print('could not plot {} with params:'.format(varname))
                print(plot_args)
                print(m)

        tab = dcc.Tab(label = model_name, children = tab_children)
        
        tabs_children.append(tab)

    tabs = dcc.Tabs(children = tabs_children)
    app.layout.children.append(tabs)
    
    return app

def generate_range_slider_callback(app, model, model_name, varname, **kwargs):
    output = Output("graph-{}-{}".format(model_name, varname), "figure")
    input_dict = dict()
    for k,v in kwargs.items():
        input_dict[k] = [
            Input("rangeslider-{}-{}-{}".format(model_name, varname, k), "value"),
            Input("input-{}-{}-{}".format(model_name, varname, k), "value")]

    inputs  = []
    for inputs_ in input_dict.values():
        inputs.extend(inputs_)

    @app.callback( output, inputs)
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


def get_range_slider(model_name, varname, **kwargs):
    divs = []
    for k,v in kwargs.items():
        marks = {i: "{0:0.2f}".format(i) for i in v[::int(len(v)/31)]}
        # print(len(marks))
        # print(marks)

        row_children = [
            html.Div(
                className = "eight columns",
                children = [
                        dcc.RangeSlider(
                            id = "rangeslider-{}-{}-{}".format(model_name, varname, k),
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
                        id = "input-{}-{}-{}".format(model_name, varname, k),
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
