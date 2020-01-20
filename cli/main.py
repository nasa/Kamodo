import hydra
import numpy as np
from plotly.offline import plot
import plotly.graph_objs as go
from omegaconf import OmegaConf
from os import path


def eval_config(params):
    """generate arrays dict from configuration"""
    args = {}
    for k,v in params.items():
        try:
            args[k] = np.array(v)
        except:
            args[k] = np.linspace(v['min'], v['max'], v['n'])
    return args

def write_plot_div(plot_result, plot_conf):
    """writes plot div to file"""
    if plot_result is not None:
        plot_filename = plot_conf['filename']
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

@hydra.main(config_path='conf/config.yaml', strict = True)
def main(cfg):
    """A low-coding command line interface for Kamodo

    This application allows users to work with kamodo-compatible models and data
    directly from the command line. 

    Custom models, data, and expressions may be composed by editing config files
    without needing to write python.
    """
    cfg = config_override(cfg)

    model = hydra.utils.instantiate(cfg.model)

    if cfg.verbose:
        print(model.detail())   
        print(model.to_latex())

    if cfg.model.evaluate is not None:
        for varname, params in cfg.model.evaluate.items():
            units = model[varname].meta['units']
            lhs = "{}({})".format(varname, ','.join(['{}={}'.format(k,v) for k,v in params.items()]))

            rhs = model[varname](**eval_config(params))

            print("{} {} = \n".format(lhs,units))
            print(rhs)

    plot_conf = cfg.model.plot_conf
    fig_layout = cfg.model.fig_layout
    plot_params = cfg.model.plot
    if plot_params is not None:
        if hasattr(plot_params, 'items'):
            pass
        elif type(plot_params) is str:
            plot_params = {plot_params: {}}
        else:
            try:
                plot_params = {k : {} for k in plot_params}
            except:
                raise IOError("cannot handle plot parameters of type {}".format(type(plot_params)))

        for varname, params in plot_params.items():
            plot_args = {varname: eval_config(params)}
            try:
                fig = go.Figure(model.plot(**plot_args))
                if fig_layout is not None:
                    fig.update_layout(**fig_layout)
                plot_result = plot(fig, **plot_conf)
                if plot_result is not None:
                    write_plot_div(plot_result, plot_conf)
            except:
                print('could not plot {} with params:'.format(varname))
                print(plot_args)
                raise


# entrypoint for package installer
def entry():
    main()

if __name__ == "__main__":
    main()