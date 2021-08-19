import hydra
import numpy as np
from plotly.offline import plot
import plotly.graph_objs as go
from omegaconf import OmegaConf
from os import path
from kamodo import Kamodo, compose
import os
from hydra import utils


def eval_config(params):
    """generate arrays dict from configuration"""
    args = {}
    if hasattr(params, 'items'):
        for k,v in params.items():
            try:
                args[k] = np.array(v)
            except:
                args[k] = np.linspace(v['min'], v['max'], v['n'])
    return args

def write_plot_div(plot_result, plot_filename):
    """writes plot div to file"""
    if plot_result is not None:
        with open(plot_filename, 'w') as f:
            f.write(plot_result)
            f.write('') # needs a newline or else embedding breaks

def config_override(cfg):
    """Overrides with user-supplied configuration

    kamodo will override its configuration using
    kamodo.yaml if it is in the current working directory
    or users can set an override config:
        config_override=path/to/myconfig.yaml
    """
    override_path = hydra.utils.to_absolute_path(cfg.config_override)
    if path.exists(override_path):
        override_conf = OmegaConf.load(override_path)
        # merge overrides first input with second
        cfg = OmegaConf.merge(cfg, override_conf)
    return cfg

def get_plot_dict(plot_params):
    if hasattr(plot_params, 'items'):
        pass
    elif type(plot_params) is str:
        plot_params = {plot_params: {}}
    else:
        try:
            plot_params = {k : {} for k in plot_params}
        except:
            raise IOError("cannot handle plot parameters of type {}".format(type(plot_params)))
    return plot_params


@hydra.main(config_path='conf/config.yaml', strict = False)
def main(cfg):
    """A low-coding command line interface for Kamodo

    This application allows users to work with kamodo-compatible models and data
    directly from the command line. 

    Custom models, data, and expressions may be composed by editing config files
    without needing to write python.
    """
    cfg = config_override(cfg)
    if cfg.verbose > 1:
        print(cfg.pretty())

    plot_conf = cfg.plot_conf
    plot_out_prefix = plot_conf.filename.split('.html')[0]
    plot_conf.pop('filename')

    models = dict()
    for model_name, model_conf in cfg.models.items():
        print(model_name)
        # model_conf['symbol_dict'] = locals()

        model = hydra.utils.instantiate(model_conf)
        models[model_name] = model

        if cfg.verbose > 0:
            print(model.detail())
        if cfg.verbose > 1:   
            print(model.to_latex())

        if model_conf.evaluate is not None:
            for varname, params in model_conf.evaluate.items():
                try:
                    units = model[varname].meta['units']
                    lhs = "{}({})".format(varname, ','.join(['{}={}'.format(k,v) for k,v in params.items()]))

                    rhs = model[varname](**eval_config(params))

                    print("{} {} = \n".format(lhs,units))
                    print(rhs)
                except Exception as m:
                    print('could not evaluate {} with params\n{}'.format(
                        varname, params))
                    print(m)
                    pass

        fig_layout = model_conf.fig_layout
        plot_params = model_conf.plot
        if plot_params is not None:
            plot_params = get_plot_dict(plot_params)

            for varname, params in plot_params.items():
                try:
                    plot_args = {varname: eval_config(params)}
                except:
                    print('issue with these params for {}'.format(varname))
                    print(params)
                    raise
                try:
                    fig = go.Figure(model.plot(**plot_args))
                    if fig_layout is not None:
                        fig.update_layout(**fig_layout)
                    out_filename = '{}_{}'.format(plot_out_prefix, varname)
                    plot_result = plot(fig, filename = out_filename, **plot_conf)
                    if plot_result is not None:
                        write_plot_div(plot_result, out_filename)
                except:
                    print('could not plot {} with params:'.format(varname))
                    print(plot_args)
                    print(plot_conf)
                    raise

    kamodo = compose(**models)
    if cfg.verbose > 0:
        print(kamodo.detail())
    if cfg.verbose > 1:
        print(kamodo.to_latex())


# entrypoint for package installer
def entry():
    main()

if __name__ == "__main__":
    main()