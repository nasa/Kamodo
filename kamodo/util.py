"""
Copyright © 2017 United States Government as represented by the Administrator, National Aeronautics and Space Administration.  
No Copyright is claimed in the United States under Title 17, U.S. Code.  All Other Rights Reserved.
"""
import logging
import os
import tempfile
import sys
import numpy.f2py  # just to check it presents
from numpy.distutils.exec_command import exec_command



from collections import OrderedDict, defaultdict
import collections
import functools
from sympy import sympify as parse_expr
from sympy.utilities.autowrap import ufuncify
import functools
from decorator import decorator, decorate
from sympy import symbols, Symbol
from sympy.core.function import UndefinedFunction
from inspect import getargspec, getfullargspec
import inspect
from sympy.physics import units
from sympy.physics import units as sympy_units
import numpy as np
from sympy import latex, Eq
from sympy.parsing.latex import parse_latex
import pandas as pd
from scipy.integrate import solve_ivp


def get_unit_quantity(name, base, scale_factor, abbrev=None, unit_system='SI'):
    '''Define a unit in terms of a base unit'''
    u = units.Quantity(name, abbrev=abbrev)
    base_unit = getattr(sympy_units, base)

    try:
        # sympy >= 1.5 raises the following warning:
        #   Use unit_system.set_quantity_dimension or
        # <unit>.set_global_relative_scale_factor
        u.set_global_relative_scale_factor(scale_factor, base_unit)
    except AttributeError:
        u.set_dimension(base_unit.dimension)
        u.set_scale_factor(scale_factor * base_unit, unit_system=unit_system)

    return u


unit_subs = dict(nT=get_unit_quantity('nanotesla', 'tesla', .000000001, 'nT', 'SI'),
                 R_E=get_unit_quantity('earth radii', 'km', 6.371e6, 'R_E', 'SI'),
                 erg=get_unit_quantity('erg', 'J', .0000001, 'erg', 'SI'),
                 nPa=get_unit_quantity('nanopascals', 'pascal', .000000001, 'nPa', 'SI'),
                 )

# global_ureg = UnitRegistry()
# global_ureg.define('m^3 = m * m * m = m3')
# global_ureg.define('cc = cm * cm * cm')
# global_ureg.define('R_s = 6.957e8 m')
# global_ureg.define('AU = 1.496e+11 m')



# def compile_fortran(source, module_name, extra_args='', folder='./'):
#     with tempfile.NamedTemporaryFile('w', suffix='.f90') as f:
#         f.write(source)
#         f.flush()

#         args = ' -c -m {} {} {}'.format(module_name, f.name, extra_args)
#         command = 'cd "{}" && "{}" -c "import numpy.f2py as f2py;f2py.main()" {}'.format(folder, sys.executable, args)
#         status, output = exec_command(command)
#         return status, output, command


def substr_replace(name, name_maps):
    """replaces all substrings in name with those given by name_maps"""
    for old, new in name_maps:
        name = name.replace(old, new)
    return name


# This should be in yaml
def beautify_latex(s):
    return substr_replace(s, [
        ('plus', '+'),
        ('minus', '-'),
        ('comma', ','),
        ('LEFT', '\\left ('),
        ('RIGHT', '\\right )'),
        ('integral', '\\int'),
        ('rvert', '\\rvert'),
        ('lvert', '\\lvert'),
    ])


def arg_to_latex(arg):
    return beautify_latex(latex(Symbol(arg)))


def decorator_wrapper(f, *args, **kwargs):
    """Wrapper needed by decorator.decorate to pass through args, kwargs"""
    return f(*args, **kwargs)


def kamodofy(_func=None, units='', data=None, update=None, equation=None, citation=None, hidden_args=[], **kwargs):
    """Adds meta and data attributes to functions for compatibility with Komodo

    meta: a dictionary containing {units: <str>}
    data:
        if supplied, set f.data = data
        if not supplied, set f.data = f(), assuming it can be called with no arguments.
            If f cannot be called with no arguments, set f.data = None
    """

    def decorator_kamodofy(f):
        f.meta = dict(units=units, citation=citation, equation=equation, hidden_args=hidden_args)
        if citation is not None:
            f.__doc__ = f.__doc__ + '\n\ncitation: {}'.format(citation)
        f.update = update
        if data is None:
            try:
                f.data = f(**kwargs)
            except:
                f.data = None
        else:
            f.data = data

        if equation is not None:
            latex_str = equation.strip("$")
            f._repr_latex_ = lambda: latex_str
        # f._repr_latex_ = lambda : "${}$".format(latex(parse_latex(latex_str)))
        else:
            f_ = symbols(f.__name__, cls=UndefinedFunction)
            lhs = f_.__call__(*symbols(getfullargspec(f).args))
            lambda_ = symbols('lambda', cls=UndefinedFunction)
            latex_eq = latex(Eq(lhs, lambda_(*lhs.args)))
            f._repr_latex_ = lambda: "${}$".format(latex(latex_eq))

        return decorate(f, decorator_wrapper)  # preserves signature

    if _func is None:
        return decorator_kamodofy
    else:
        return decorator_kamodofy(_func)


def symbolic(sym_names):
    """Returns symbolic functions"""
    syms = symbols(sym_names, cls=UndefinedFunction)
    if len(syms) == 1:
        return syms[0]
    else:
        return syms


def sort_symbols(symbols):
    symbols_ = list(symbols)
    symbols_.sort(key=str)
    return tuple(symbols_)


def valid_args(f, kwargs):
    '''Extract arguments from kwargs that appear in f'''
    valid = OrderedDict()
    if type(f) == np.vectorize:
        for a in getargspec(f.pyfunc).args:
            if a in kwargs:
                valid[a] = kwargs[a]
        return valid
    else:
        for a in getfullargspec(f).args:
            if a in kwargs:
                valid[a] = kwargs[a]
        return valid


def eval_func(func, kwargs):
    '''Evaluate function over valid argument'''
    try:
        return func(**valid_args(func, kwargs))
    except TypeError as m:
        raise TypeError(str(m) + str(list(valid_args(func, kwargs).keys())))


def get_defaults(func):
    sig = inspect.signature(func)
    defaults = {}
    for k, v in sig.parameters.items():
        if v.default is not inspect._empty:
            defaults[k] = v.default

    return defaults

    spec = getargspec(func)
    args = spec.args
    defaults = spec.defaults
    if defaults is not None:
        return OrderedDict(list(zip(args[-len(defaults):], spec.defaults)))
    else:
        return None


def cast_0_dim(a, to):
    if a.ndim == 0:
        return a * np.ones(to.shape)
    else:
        return a


def simulate(state_funcs, **kwargs):
    """Iterate over functions.

    state_funcs(OrderedDict)
        key: the variable to update
        value: the function that updates the variable

    Any remaining kwargs are passed to the state functions.

    The order of the keys in state_funcs determines which variables get updated first.
    """

    verbose = kwargs.get('verbose', False)
    steps = kwargs.get('steps', 1)

    state_dict = kwargs.copy()
    yield OrderedDict([(k, state_dict.get(k, None)) for k in state_funcs])
    for i in range(steps):
        result = []
        for arg, f in list(state_funcs.items()):
            try:
                result.append((arg, eval_func(f, state_dict)))
            except TypeError as m:
                raise TypeError('{}:'.format(arg) + str(m))
            state_dict.update(OrderedDict(result))
        yield OrderedDict(result)


def pad_nan(array):
    if len(array.shape) == 1:
        # print('padding array {}'.format(array.shape))
        return np.pad(array.astype(float), (0, 1),
                      mode='constant', constant_values=np.nan).T
    elif len(array.shape) == 2:
        # print('padding array {}'.format(array.shape))
        return np.pad(array.astype(float), ((0, 1), (0, 0)),
                      mode='constant', constant_values=np.nan)
    else:
        raise NotImplementedError('cannot pad shape {}'.format(array.shape))


def concat_solution(gen, variable):
    result = []
    params = defaultdict(list)
    for f in gen:
        for k, v in list(get_defaults(f).items()):
            params[k].append(pad_nan(v))
        params[variable].append(pad_nan(f.data))

    for k, v in list(params.items()):
        if len(v[0].shape) == 1:
            params[k] = np.hstack(v)
        else:
            params[k] = np.vstack(v)
    return params


existing_plot_types = pd.DataFrame({
    ('1', 'N', 'N'): ['1-d line', ''],
    ('1', 'N', 'Nx2'): ['2-d line', ''],
    ('1', 'N', 'Nx3'): ['3-d line', ''],
    ('3', 'N, N, N', 'N'): ['3-d colored line', ''],
    ('1', 'Nx2', 'Nx2'): ['2-d vector field', ''],
    ('1', 'Nx3', 'Nx3'): ['3-d vector field', ''],
    ('2', 'N, M', 'NxM'): ['2-d contour', 'indexing'],
    ('2', 'NxM, NxM', 'NxM'): ['2-d contour (skew/carpet)', 'indexing'],
    ('3', 'NxM, NxM, NxM', '1'): ['Parametric Surface', ''],
    ('3', 'NxM, NxM, NxM', 'NxM'): ['Coloured Parametric Surface', ''],
    ('3', 'L, M, 1', 'LxM'): ['Map-to-plane', 'indexing*'],
    ('3', '1, M, N', 'MxN'): ['Map-to-plane', 'indexing*'],
    ('3', 'L, 1, N', 'LxN'): ['Map-to-plane', 'indexing*']
}).T
existing_plot_types.index.set_names(['nargs', 'arg shapes', 'out shape'], inplace=True)
existing_plot_types.columns = ['Plot Type', 'notes']

# manually generate the appropriate function signature
grid_wrapper_def = r"""def wrapped({signature}):
    coordinates = np.meshgrid({arg_str}, indexing = 'xy')
    points = np.column_stack([c.ravel() for c in coordinates])
    return np.squeeze({fname}(points).reshape(coordinates[0].shape, order = 'A'))
    """


def gridify(_func=None, **defaults):
    """Given a function of shape (n,dim) and arguments of shape (L), (M), calls f with points L*M"""

    def decorator_gridify(f):

        arg_str = ', '.join([k for k in defaults])

        signature = ''
        for k, v in defaults.items():
            signature = signature + "{} = {},".format(k, k)

        scope = {**defaults}
        scope['np'] = np
        scope[f.__name__] = f

        exec(grid_wrapper_def.format(signature=signature, arg_str=arg_str, fname=f.__name__), scope)
        wrapped = scope['wrapped']
        wrapped.__name__ = f.__name__
        wrapped.__doc__ = f.__doc__

        return decorate(wrapped, decorator_wrapper)

    if _func is None:
        return decorator_gridify
    else:
        return decorator_gridify(_func)


def pointlike(_func=None, signature=None, otypes=[np.float], squeeze=None):
    """Transforms a single-argument function to one that accepts m points of dimension n"""

    def decorator_pointlike(func):
        def argument_wrapper(f, *args, **kwargs):
            """Wrapper needed by decorator.decorate to pass through args, kwargs"""
            if type(args[0]) != np.array:
                args = [np.array(x) for x in args]
            for i, x in enumerate(args):
                if len(x.shape) == 1:
                    args[i] = np.expand_dims(x, axis=0)
            if squeeze is not None:
                try:
                    return np.vectorize(f, otypes=otypes, signature=signature)(*args, **kwargs).squeeze(squeeze)
                except:
                    return np.vectorize(f, otypes=otypes, signature=signature)(*args, **kwargs)
            else:
                return np.vectorize(f, otypes=otypes, signature=signature)(*args, **kwargs)

        if not hasattr(func, '__name__'):
            func.__name__ = 'pointlike'

        return decorate(func, argument_wrapper)

    if _func is None:
        return decorator_pointlike
    else:
        return decorator_pointlike(_func)


def event(func, terminal=True, direction=0):
    def wrapped(*args, **kwargs):
        return func(*args, **kwargs)

    wrapped.terminal = terminal
    wrapped.direction = direction
    return wrapped


def solve(fprime=None, seeds=None, varname=None, interval=None,
          dense_output=True,  # generate a callable solution
          events=None,  # stop when event is triggered
          vectorized=True,
          npoints=50,
          directions=(-1, 1),
          verbose=False,
          ):
    """Decorator that solves initial value problem for a given function

    Can be used to generate streamlines, streaklines, fieldlines, etc
    """

    if len(seeds.shape) > 1:
        pass
    else:
        seeds = np.expand_dims(seeds, axis=0)

    t_eval = np.linspace(*interval, npoints)
    nseeds = len(seeds)

    def decorator_solve(f):
        solutions = []
        t = []

        fprime_ = {}
        for d in directions:
            fprime_[d] = lambda s, y: d * f(y.T)

        for i, seed in enumerate(seeds):
            for d in directions:
                result = solve_ivp(fprime_[d], interval, seed,
                                   dense_output=dense_output,
                                   events=events,
                                   vectorized=vectorized,
                                   t_eval=t_eval)
                solutions.append(result['sol'])
                interval_bounded = result['t']
                seed_numbers = np.ones(len(interval_bounded)) * i  # *len(directions) + 1*(d > 0)
                integrals = interval_bounded[::d] * d * 1j
                if d < 0:
                    t.extend(list(seed_numbers + integrals))
                else:
                    t.extend(list(seed_numbers + integrals)[1:])

        t = np.hstack(t)

        def solution(s=t):
            s = np.array(s)
            if len(s.shape) == 0:
                s = np.expand_dims(s, axis=0)

            isolution = np.floor(s.real).astype(int) * len(directions) + (s.imag > 0)

            results = []
            seed_number = []
            integral = []
            for soln, imag_ in zip(isolution, s.imag):
                seed_number.append(np.floor(soln / len(directions)))
                integral.append(imag_)
                try:
                    if np.isnan(abs(soln + imag_)):
                        results.append(np.ones(isolution.shape[-1]) * np.nan)
                    else:
                        results.append(solutions[soln](np.abs(imag_)))
                except:
                    results.append(np.ones(isolution.shape[-1]) * np.nan)
            index_ = pd.MultiIndex.from_arrays([seed_number, integral],
                                               names=['seed', 'integral'])
            return pd.DataFrame(np.vstack(results), index=index_).drop_duplicates()

        solution.__name__ = varname

        return decorate(solution, decorator_wrapper)  # preserves signature

    if fprime is None:
        return decorator_solve
    else:
        return decorator_solve(fprime)

