import hydra
from hydra.experimental import compose
from omegaconf import OmegaConf
import dash
import dash_core_components as dcc
import dash_html_components as html
from kamodo.cli.main import eval_config
from kamodo import get_defaults
import plotly.graph_objs as go
from dash.dependencies import Input, Output, ClientsideFunction, State
import numpy as np
from os import path
import os
import numpy as np
from sympy.core.function import UndefinedFunction 
import pandas as pd
import sys 
from dash_katex import DashKatex
from plotly.subplots import make_subplots
from dash.exceptions import PreventUpdate
from flask_restful import reqparse, abort, Api, Resource
import flask
import json
from kamodo import getfullargspec
from flask.json import jsonify
from flask import request
from io import StringIO


try:
    hydra.experimental.initialize(strict=False)
except:
    pass

from ast import literal_eval


def get_equation_divs(model, model_name, model_params):

    equations = []

    options = []
    selected = []
    for var_symbol, func in model.items():
        if type(var_symbol) != UndefinedFunction:
            var_label = str(type(var_symbol))

            # make sure user-specified params get a check mark
            if str(type(var_symbol)) in model_params:
                selected.append(var_label)

            var_docs = model[var_symbol].__doc__.split('\n')
            doc_summary = var_docs[0]
            doc_complete = ''
            if len(var_docs) > 1:
                doc_complete = [html.Pre(doc_) for doc_ in var_docs[1:]]
            options.append({'label': '', 'value': var_label})
            equations.append(
                html.Div([
                    html.Div(
                        DashKatex(
                            id = "{}-{}-expression".format(model_name, str(var_symbol)),
                            expression = "{}".format(model.to_latex(
                                keys = [str(var_symbol)],
                                mode = "plain")),
                            ),
                        className = 'four columns'),
                    
                    html.Details([
                        html.Summary(doc_summary),
                        html.Div(children=doc_complete)],
                        className = 'eight columns',
                        style = {'text-align': 'left'})
                    ],
                    className = 'row'
                    )
                )

    equation_divs = html.Div([
            dcc.Checklist(
                id = 'checklist-{}'.format(model_name),
                options = options,
                value = selected,
                className = "one columns",
                style = {'text-align': 'right'},
                ),
            html.Div(
                children = equations,
                className = "eleven columns"
                ),
            ],
        className = 'row',
        )

    return equation_divs

def get_gui(cfg, server):

    app = dash.Dash(__name__,
        server=server,
        routes_pathname_prefix = cfg.gui.routes_pathname_prefix,
        external_stylesheets = cfg.gui.external_stylesheets,
        external_scripts = cfg.gui.external_scripts,
        )
    app.title = "kamodo"
    app.config.suppress_callback_exceptions = False

    tabs_children, graphs = get_tabs_children(cfg, app)

    app.layout = html.Div(children=[
        html.Div(
            [
            html.Div(
                html.Img(
                src="assets/Kamodo-black.png",
                alt="Kamodo Software Suite",
                style={'height':'80px'},
                ),
                className = 'six columns'
                ),
            html.Div(
                html.Img(
                src="assets/ccmc_logo.svg",
                alt="Kamodo Software Suite",
                style={'height':'80px'},
                ),
                className = 'six columns',
                style = {'text-align': 'right'}, 
                ),
            ],
            className = 'row'
            ),
        

        html.Div(children='''
            A low-coding interface for models and data
        '''),
        html.Div(children = [
            dcc.Textarea(
                id = "kamodo-config",
                title = "Configuration for Kamodo",
                placeholder = 'Supply a configuration for Kamodo',
                value = "{}".format(cfg.pretty()),
                style = {'width': '500px', 'height':'400px'},
                className = 'three columns'),
            html.Div(
                id = 'config-render',
                className = 'three columns',
                style = {'height': '400px'}
                ),
            dcc.Input(
                id = 'save-as',
                value = cfg.config_override,
                placeholder = "enter filename",
                className = 'two columns'
                ),
            html.Button(
                'Save',
                id = 'save-button',
                n_clicks = 0,
                )
            ],
            className = 'row',
            style = {'overflow': 'auto'}
        ),
        dcc.Tabs(id = 'kamodo-tabs',
            children = tabs_children),
    ])


    generate_config_render_callback(app)

    generate_kamodo_tabs_callback(app)

    generate_save_button_callbacks(app)

    for graph_id, (model, model_params, model_name) in graphs.items():
        generate_subplot_callback(app, graph_id, model, model_name, model_params)


    
    return app

def generate_save_button_callbacks(app):
    @app.callback(
        Output('save-button', 'children'),
        [Input('save-button', 'n_clicks'), ],
        [State('save-as', 'value'), State('kamodo-config', 'value')]
        )
    def save_config(n_clicks, filename, conf):
        if n_clicks > 0:
            try:
                cfg = OmegaConf.create(conf)
                with open(filename, 'w') as f:
                    f.write(cfg.pretty())
                    return 'saved'
            except:
                return html.Div('error', style = {'color':'red'})
        else:
            raise PreventUpdate
            return 'save'

    @app.callback(
        Output('save-as','value'),
        [Input('kamodo-config', 'value')])
    def update_save_value(conf):
        try:
            cfg = OmegaConf.create(conf)
            filename = cfg.config_override
            return filename
        except:
            raise PreventUpdate

    @app.callback(
        Output('save-button', 'n_clicks'),
        [Input('kamodo-config', 'value')])
    def update_save_button(conf):
        try:
            OmegaConf.create(conf)
            return 0
        except:
            raise PreventUpdate



def generate_config_render_callback(app):
    @app.callback(
        Output('config-render', 'children'),
        [Input('kamodo-config', 'value')])
    def config_render(conf):
        try:
            cfg = OmegaConf.create(conf)
            return dcc.Markdown(
                id = 'config-render-markdown',
                children = "```yaml\n{}\n```".format(cfg.pretty()),
                )
        except Exception as m:
            return html.Div(str(m), style = {'color':'red'})



def get_tabs_children(cfg, app):

    tabs_children = []
    graphs = {}

    for model_name, model_conf in cfg.models.items():
        print(model_name)
        tab_children = []
        tab_children.append(html.Br())

        try:
            model = hydra.utils.instantiate(model_conf)
        except Exception as m:
            print('could not instantiate {}'.format(model_name))
            raise

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
     

        tab_children.append(get_equation_divs(model, model_name, model_params))

        graph_id = "graph-{}".format(model_name)

        tab_children.append(
            html.Div(
                html.Div(
                    dcc.Graph(
                        id = graph_id,
                        # figure = fig_subplots,
                        style = {
                            'height': '1000px', 
                            'width':'1200px', 
                            'text-align': 'center',
                            'margin': 'auto'}),
                    className = 'twelve columns'),
                className = 'row'),
            )

        try:
            plots = get_kamodo_plots(model, model_params.keys(), model_params)
        except Exception as m:

            tab_children.append(html.Div(str(m), style = {'color': 'red'}))
            plots = {}
                    
        tab_children.append(
            dcc.Store(id = graph_id + '-store',
                data = plots),
            )



        # fig_subplots = make_kamodo_subplots(model, model_params.keys(), model_params)
 



        graphs[graph_id] = model, model_params, model_name

        tab = dcc.Tab(label = model_name, children = tab_children)
        
        tabs_children.append(tab)

    return tabs_children, graphs

def generate_kamodo_tabs_callback(app):
    @app.callback(Output('kamodo-tabs','children'),
    [Input('kamodo-config', 'value')])
    def kamodo_content(conf):
        try:
            cfg = OmegaConf.create(conf)
            tabs_children, graphs = get_tabs_children(cfg, app)
            return tabs_children
        except Exception as m:
            print("could not get tabs children")
            print(m)
            print(conf)
            return html.Pre(str(m).split('\n'), style = {'color':'red'})
            

def get_kamodo_plots(model, var_names, model_params):
    plots = {}
    for var_symbol, func in model.items():
        if type(var_symbol) != UndefinedFunction:
            var_label = str(type(var_symbol))

            fig = model.figure(
                var_label,
                indexing = 'ij',
                return_type = True, 
                **model_params.get(var_label, {}))
            plots[var_label] = fig['data'][0]

    return plots


# def make_kamodo_subplots(model, var_names, model_params):
#     fig_subplots = make_subplots(rows=len(var_names), cols=1)
#     for plot_count, varname in enumerate(var_names):
#         try:
#             fig = model.figure(
#                 varname,
#                 indexing = 'ij',
#                 return_type = True, 
#                 **model_params.get(varname, {}))
#             fig_subplots.add_trace(fig['data'][0], row = plot_count+1, col = 1)
#         except:
#             print("could not plot {} with model_params:".format(varname))
#             print(model_params)
#             raise PreventUpdate
#     return fig_subplots


def generate_subplot_callback(app, graph_id, model, model_name, model_params):
    app.clientside_callback(
        ClientsideFunction("clientside", "subplots"),
        Output(component_id=graph_id, component_property="figure"),
        [Input(graph_id + '-store', "data"),
        Input('checklist-{}'.format(model_name), 'value')],
    )


    

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

def get_api(app):
    api = Api(app)
    return api


def main():

    cfg = compose('conf/kamodo.yaml')

    cli_conf = OmegaConf.from_cli()
    
    cfg = OmegaConf.merge(cfg, cli_conf)
        
    config_override = None
    extra_files = []
    if cfg.config_override is not None:
        override_path = "{}/{}".format(os.getcwd(),cfg.config_override)
        if path.exists(override_path):
            print("found {}".format(override_path))
            config_override = OmegaConf.load(override_path)
            # print(config_override.pretty())
            extra_files.append(override_path)
        else:
            if cfg.verbose > 0:
                print("could not get override: {}".format(override_path))
    

        if config_override is not None:
            cfg = OmegaConf.merge(cfg, config_override)
    

    if cfg.verbose > 0:
        print(cfg.pretty())

    server = flask.Flask(__name__)

    # app = Flask(__name__)
    api = Api(server)

    models = dict()
    for model_name, model_conf in cfg.models.items():
        try:
            # this might fail
            model_ = hydra.utils.instantiate(model_conf)
            models[model_name] = model_
        except:
            pass
    
    class Models(Resource):
        def get(self):
            details = dict()
            for model_name, model_ in models.items():
                detail = model_.detail().astype(str)
                print(detail)
                details[model_name] = detail.to_dict(
                    # default_handler=str, 
                    # indent =4, 
                    orient = 'index',
                    )
            return details

    api.add_resource(Models, '/api', '/api/')


    for model_name, model_ in models.items():
        api.add_resource(
            get_model_resource(model_name, model_),
            '/api/{}'.format(model_name),
            '/api/{}/'.format(model_name),
            endpoint = model_name)

        for var_symbol in model_:
            if type(var_symbol) != UndefinedFunction:
                var_label = str(type(var_symbol))
                api.add_resource(
                    get_func_resource(model_name, model_, var_symbol),
                    '/api/{}/{}'.format(model_name, var_label),
                    endpoint =  '/'.join([model_name, var_label]))



    @server.route('/')
    def index():
        return 'Hello Flask app'

    app = get_gui(cfg, server)
    app.run_server(debug = True, extra_files = extra_files)


def get_model_resource(model_name, model):
    class modelResource(Resource):
        def get(self):
            return model.detail().astype(str).to_dict(orient='index')

    return modelResource

def get_func_resource(model_name, model, var_symbol):
    parser = reqparse.RequestParser()
    func = model[var_symbol]
    defaults = get_defaults(func)
    for arg in getfullargspec(func).args:
        if arg in defaults:
            parser.add_argument(arg, type = str)
        else:
            parser.add_argument(arg, type = str, required = True)

    class funcResource(Resource):
        def get(self):
            
            args_ = parser.parse_args(strict = True)
            args = dict()

            for argname, val_ in args_.items():
                args[argname] = pd.read_json(StringIO(val_), typ = 'series')

            return func(**args).tolist()
      
            

    return funcResource

# entrypoint for package installer
def entry():
    main()


if __name__ == "__main__":
    main()
