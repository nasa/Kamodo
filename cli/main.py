import hydra
import numpy as np
from plotly.offline import plot
import plotly.graph_objs as go


def eval_config(params):
    args = {}
    for k,v in params.items():
        try:
            args[k] = np.array(v)
        except:
            args[k] = np.linspace(v['min'], v['max'], v['n'])
    return args


@hydra.main(config_path='conf/config.yaml', strict = False)
def main(cfg):
    model = hydra.utils.instantiate(cfg.model)
    print(model.detail())
    print()
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
            plot(fig, **plot_conf)
        except:
            print('could not plot {} with params:'.format(varname))
            print(plot_args)


def entry():
    main()

if __name__ == "__main__":
    main()