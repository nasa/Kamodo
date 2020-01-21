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

hydra.experimental.initialize()


def get_gui(cfg):
    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

    external_scripts = ['https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.4/MathJax.js?config=TeX-MML-AM_CHTML']


    app = dash.Dash(__name__, 
        external_stylesheets = external_stylesheets,
        external_scripts = external_scripts)


    app.layout = html.Div(children=[
        html.H1(children='Kamodo'),

        html.Div(children='''
            Kamodo: A low-coding interface for models and data.
        '''),
    ])

    for model_name, model_conf in cfg.models.items():
        print(model_name)

        try:
            model = hydra.utils.instantiate(model_conf)
        except Exception as m:
            print('could not instantiate {}'.format(model_name))
            print(m)


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
                # print(getattr(fig,'chart_type', 'no chart type'))
                app.layout.children.append(dcc.Graph(
                    id = "graph-{}-{}".format(model_name, varname),
                    figure = go.Figure(fig),
                    ))
                if chart_type == 'line':
                    app.layout.children.extend(get_range_slider(
                        model_name, varname, **plot_args[varname]))
                    generate_range_slider_callback(
                        app,
                        model,
                        model_name,
                        varname,
                        **plot_args[varname])

            except Exception as m:
                print('could not plot {} with params:'.format(varname))
                print(plot_args)
                print(m)


    
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
        divs.append(
            dcc.RangeSlider(
                id = "rangeslider-{}-{}-{}".format(model_name, varname, k),
                marks=marks,
                min=min(v),
                max=max(v),
                step = (max(v) - min(v))/len(v),
                value=[min(v), max(v)]
            ) )
        divs.append(
            dcc.Input(
                id = "input-{}-{}-{}".format(model_name, varname, k),
                placeholder='number of points in {}'.format(k),
                type='number',
                value = len(v),
                min = 3,
                max = 10*len(v),
            )  )
    return divs

def main():

    cfg = compose('conf/kamodo.yaml', overrides=["config_override=../../kamodo.yaml"])

    # print(cfg.pretty())
    cfg_overrides = compose(cfg.config_override)
    # print(cfg_overrides.pretty())    
    cfg = OmegaConf.merge(cfg, cfg_overrides)
    # print(cfg.pretty())
    print(cfg.plot_conf.pretty())

    app = get_gui(cfg)
    app.run_server(debug = True)

# entrypoint for package installer
def entry():
    main()


if __name__ == "__main__":
    main()
