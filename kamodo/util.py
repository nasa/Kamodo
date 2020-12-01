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
from sympy.physics.units.util import _get_conversion_matrix_for_expr

from sympy.core.compatibility import reduce, Iterable, ordered
from sympy import Add, Mul, Pow, Tuple, sympify, default_sort_key
from sympy.physics.units.quantities import Quantity

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


unit_subs = dict(R_E=get_unit_quantity('earth radii', 'km', 6.371e6, 'R_E', 'SI'),
                 erg=get_unit_quantity('erg', 'J', .0000001, 'erg', 'SI'),
                 )

sympy_units.erg = unit_subs['erg']

prefix_dict = sympy_units.prefixes.PREFIXES  #built-in dictionary of Prefix instances
#test_unit_subs={}  #dictionary to replace subs in kamodo.get_unit_quantities()
unit_list = ['m', 's', 'g', 'A', 'K', 'radian', 'sr', 'cd', 'mole', 'eV', 'Pa', 'F', 'N',
             'V', 'Hz', 'C', 'W', 'Wb', 'H', 'S', 'Bq', 'Gy', 'erg', 'T']

#list of SI units included in sympy (likely not complete)

for item in unit_list:
    unit_item = getattr(sympy_units, item)
    unit_subs[item] = unit_item
    for key in prefix_dict.keys():
        unit_subs[key+item] = get_unit_quantity(str(prefix_dict[key].name)+
                                                str(unit_item.name), str(unit_item.name),
                                                float(prefix_dict[key].scale_factor),
                                                abbrev=key+item)



def substr_replace(name, name_maps):
    """replaces all substrings in name with those given by name_maps"""
    for old, new in name_maps:
        name = name.replace(old, new)
    return name


# This should be in yaml
def beautify_latex(s):
    return substr_replace(s, [
        ('**', '^'),
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


def kamodofy(_func=None, units='', arg_units = None, data=None, update=None, equation=None, citation=None, hidden_args=[], **kwargs):
    """Adds meta and data attributes to functions for compatibility with Komodo

    meta: a dictionary containing {units: <str>}
    data:
        if supplied, set f.data = data
        if not supplied, set f.data = f(), assuming it can be called with no arguments.
            If f cannot be called with no arguments, set f.data = None
    """

    def decorator_kamodofy(f):
        f.meta = dict(units=units, arg_units = arg_units, citation=citation, equation=equation, hidden_args=hidden_args)
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
    coordinates = np.meshgrid({arg_str}, indexing = 'xy', sparse = False, copy = False)
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



def convert_to(expr, target_units, unit_system="SI", raise_errors=True):
    """
    Same as sympy.convert_to but accepts equations and allows functions of units to pass

    Convert ``expr`` to the same expression with all of its units and quantities
    represented as factors of ``target_units``, whenever the dimension is compatible.

    ``target_units`` may be a single unit/quantity, or a collection of
    units/quantities.

    Examples
    ========

    >>> from sympy.physics.units import speed_of_light, meter, gram, second, day
    >>> from sympy.physics.units import mile, newton, kilogram, atomic_mass_constant
    >>> from sympy.physics.units import kilometer, centimeter
    >>> from sympy.physics.units import gravitational_constant, hbar
    >>> from sympy.physics.units import convert_to
    >>> convert_to(mile, kilometer)
    25146*kilometer/15625
    >>> convert_to(mile, kilometer).n()
    1.609344*kilometer
    >>> convert_to(speed_of_light, meter/second)
    299792458*meter/second
    >>> convert_to(day, second)
    86400*second
    >>> 3*newton
    3*newton
    >>> convert_to(3*newton, kilogram*meter/second**2)
    3*kilogram*meter/second**2
    >>> convert_to(atomic_mass_constant, gram)
    1.660539060e-24*gram

    Conversion to multiple units:

    >>> convert_to(speed_of_light, [meter, second])
    299792458*meter/second
    >>> convert_to(3*newton, [centimeter, gram, second])
    300000*centimeter*gram/second**2

    Conversion to Planck units:

    >>> from sympy.physics.units import gravitational_constant, hbar
    >>> convert_to(atomic_mass_constant, [gravitational_constant, speed_of_light, hbar]).n()
    7.62963085040767e-20*gravitational_constant**(-0.5)*hbar**0.5*speed_of_light**0.5

    """


    from sympy.physics.units import UnitSystem
    unit_system = UnitSystem.get_unit_system(unit_system)

    if not isinstance(target_units, (Iterable, Tuple)):
        target_units = [target_units]


    if hasattr(expr, 'rhs'):
        return Eq(convert_to(expr.lhs, target_units, unit_system),
                 convert_to(expr.rhs, target_units, unit_system))
    if type(type(expr)) is UndefinedFunction:
        print('undefined input expr:{}'.format(expr))
        return expr

    if isinstance(expr, Add):
        return Add.fromiter(convert_to(i, target_units, unit_system) for i in expr.args)

    expr = sympify(expr)

    if not isinstance(expr, Quantity) and expr.has(Quantity):
        expr = expr.replace(lambda x: isinstance(x, Quantity), lambda x: x.convert_to(target_units, unit_system))

    def get_total_scale_factor(expr):
        if isinstance(expr, Mul):
            return reduce(lambda x, y: x * y, [get_total_scale_factor(i) for i in expr.args])
        elif isinstance(expr, Pow):
            return get_total_scale_factor(expr.base) ** expr.exp
        elif isinstance(expr, Quantity):
            return unit_system.get_quantity_scale_factor(expr)
        return expr

    depmat = _get_conversion_matrix_for_expr(expr, target_units, unit_system)
    if depmat is None:
        if raise_errors:
            raise NameError('cannot convert {} to {} {}'.format(expr, target_units, unit_system))
        return expr

    expr_scale_factor = get_total_scale_factor(expr)
    return expr_scale_factor * Mul.fromiter((1/get_total_scale_factor(u) * u) ** p for u, p in zip(target_units, depmat))

def resolve_unit(expr, unit_registry, verbose=False):
    """get the registered unit for the expression

    unit_registry {f(x): f(cm), f(cm): kg/m^2}
    """
    unit = unit_registry.get(expr, None)
    if verbose:
        print('resolve_unit: {}'.format(unit))
    for k, k_unit in unit_registry.items():
        if str(expr) == str(k):
            unit = k_unit
            continue
    if verbose:
        print('resolve_unit: after registry {}'.format(unit))

    if unit in unit_registry:
        return unit_registry[unit]

    if isinstance(unit, UndefinedFunction):
        # {f(x):g(x)}
        result = unit_registry.get(unit)
        if verbose:
            print('resolve_unit: returning {}'.format(result))
        return result
    else:
        if verbose:
            print('{} is not an UndefinedFunction'.format(unit))

    if unit is None:
        unit = expr.subs(unit_registry)

    if len(unit.free_symbols) > 0:
        return None

    if unit is not None:
        if hasattr(unit, 'dimension'):
            return unit
        if isinstance(unit, Mul):
            return unit
        if isinstance(unit, Pow):
            return unit
        return unit_registry.get(unit, None)

def get_expr_unit(expr, unit_registry, verbose=False):
    '''Get units from an expression'''

    for func in unit_registry:
        if type(expr) == type(func):
            # b(a) = b(x)
            if verbose:
                print('get_expr_unit: found matching {}:'.format(func))
                print('get_expr_unit: func free symbols {}'.format(func.free_symbols))
                print('get_expr_unit: expr free_symbols: {}'.format(expr.free_symbols))
            func_units = resolve_unit(func, unit_registry, verbose)
            arg_units = get_arg_units(func, unit_registry)
            if verbose:
                print('get_expr_unit:   matching func units:', func_units)
                print('get_expr_unit:   matching arg units:', arg_units)
                print('get_expr_unit:   expression args:', expr.args)
                print('get_expr_unit:   expression arg units:',
                      [resolve_unit(arg, unit_registry, verbose) for arg in expr.args])
            assert len(arg_units) == len(expr.args)
            arg_swap = dict()
            for func_arg, expr_arg in zip(func.args, expr.args):
                if verbose:
                    print('get_expr_unit:   from {}:{} to {}:{}'.format(
                        func_arg,
                        arg_units[func_arg],
                        expr_arg,
                        convert_to(
                            resolve_unit(expr_arg, unit_registry, verbose),
                            arg_units[func_arg])))
                arg_swap[func_arg] = expr_arg*convert_to(
                    resolve_unit(expr_arg, unit_registry, verbose),
                    arg_units[func_arg])/arg_units[func_arg]
            result = func.subs(arg_swap)
            unit_registry[result] = func_units
            return result
        else:
            if verbose:
                print('get_expr_unit: {} not a match for {}'.format(func, expr))

    expr_unit = expr.subs(unit_registry, simultaneous=False)

    if isinstance(expr_unit, Add):
        # use the first term
        arg_0 = expr_unit.args[0]
        convert_to(expr_unit, arg_0)
        result = arg_0
    else:
        result = expr_unit

    if len(result.free_symbols) > 0:
        return None

    return result

def is_undefined(expr):
    return type(type(expr)) is UndefinedFunction


def get_arg_units(expr, unit_registry):
    """for each argument, retrieve the corresponding units from registered function"""
    arg_units = dict()
    if hasattr(expr, 'args') & (expr in unit_registry):
        # f(x,y,z) in unit_registry
        unit_signature = unit_registry[expr]
        # f(cm,cm,cm)
        if unit_signature is not None:
            if is_undefined(unit_signature):
                for arg_, arg_unit in zip(expr.args, unit_registry[expr].args):
                    arg_units[arg_] = arg_unit
    return arg_units

def replace_args(expr, from_map, to_map):
    func_symbol = type(expr)
    arg_map = dict()
    for arg in expr.args:
        if (arg in from_map) & (arg in to_map):
            from_unit = from_map[arg]
            to_unit = to_map[arg]
            try:
                arg_map[arg] = convert_to(arg*to_unit, from_unit)/from_unit
            except:
                print('cannot convert', arg*to_unit, from_unit)
                print(arg, type(arg), to_unit, type(to_unit), from_unit, type(from_unit))
                raise
    return expr.subs(arg_map)

def unify(expr, unit_registry, to_symbol=None, verbose=False):
    """adds unit conversion factors to composed functions"""
    if verbose:
        print('unify: to_symbol:', to_symbol)
    if hasattr(expr, 'rhs'):
        return Eq(expr.lhs, unify(
            expr.rhs,
            unit_registry,
            to_symbol=expr.lhs,
            verbose=verbose))
        # everything below skipped
        if verbose:
            print('unify: {} has rhs'.format(expr))
        lhs_unit = resolve_unit(expr.lhs, unit_registry, verbose)
        if lhs_unit is not None:
            if verbose:
                print('unify: found lhs_unit:', lhs_unit)
            return Eq(expr.lhs, unify(expr.rhs, unit_registry, expr.lhs, verbose=verbose))
        else:
            if verbose:
                print('unify: could not find lhs_unit from {}'.format(expr.lhs))
                print('unify: ', unit_registry.keys())
        return Eq(expr.lhs, unify(expr.rhs, unit_registry, expr.lhs, verbose=verbose))

    if isinstance(expr, Add):
        if verbose:
            print('unify: Adding expression: {} -> {}'.format(
                expr, resolve_unit(to_symbol, unit_registry, verbose)))
        return Add.fromiter([unify(arg, unit_registry, to_symbol, verbose=verbose) for arg in expr.args])

    # if isinstance(expr, Mul):
    #     if verbose:
    #         print('unify: Multiplying expression: {} -> {}'.format(
    #             expr, resolve_unit(to_symbol, unit_registry)))

    expr_unit = resolve_unit(expr, unit_registry, verbose)

    if verbose:
        print('unify: expr unit {}'.format(expr_unit))

    if is_undefined(expr):
        if verbose:
            print('unify: undefined expression: {}'.format(expr))
        if to_symbol is not None:
            if verbose:
                print('unify: to_symbol args: {}'.format(to_symbol.args))
                print('unify: to_symbol free symbols: {}'.format(to_symbol.free_symbols))
                print('unify: expr args: {}'.format(expr.args))
                print('unify: expr free symbols: {}'.format(expr.free_symbols))
            for arg in set(to_symbol.free_symbols).intersection(expr.free_symbols):
                if verbose:
                    print('unify: replacing {} in {}'.format(arg, expr))
                expr_units = get_arg_units(expr, unit_registry)
                to_units = get_arg_units(to_symbol, unit_registry)
                if verbose:
                    print('unify: expression arg units', expr_units)
                    print('unify: to arg units', to_units)
                expr = replace_args(expr, expr_units, to_units)
            if verbose:
                print('unify: converted expression:', expr)

    if (to_symbol is not None) & (expr_unit is not None):
        to_unit = resolve_unit(to_symbol, unit_registry, verbose)
        if verbose:
            print('unify: to_unit {}'.format(to_unit))
        if get_dimensions(expr_unit) == get_dimensions(to_unit):
            if verbose:
                print('unify: {} [{}] -> to_symbol: {}[{}]'.format(
                    expr, expr_unit, to_symbol, to_unit))
            expr = convert_to(expr*expr_unit, to_unit)/to_unit
        else:
            if verbose:
                print('unify: registry:')
                for k, v in unit_registry.items():
                    print('unify:\t{} -> {}'.format(k, v))
            raise NameError('cannot convert {} [{}] to {}[{}]'.format(expr, expr_unit, to_symbol, to_unit))

    return expr

def get_abbrev(unit):
    """get the abbreviation for a mixed unit"""
    if hasattr(unit, 'abbrev'):
        return unit.abbrev
    if isinstance(unit, Mul):
        return Mul.fromiter([get_abbrev(arg) for arg in unit.args])
    if isinstance(unit, Pow):
        base, exp = unit.as_base_exp()
        return Pow(get_abbrev(base), exp)
    return unit


def get_dimensions(unit):
    """get the set of basis units"""
    if hasattr(unit, 'dimension'):
        return unit.dimension
    if isinstance(unit, Mul):
        return Mul.fromiter([get_dimensions(arg) for arg in unit.args])
    if isinstance(unit, Pow):
        base, exp = unit.as_base_exp()
        return Pow(get_dimensions(base), exp)
    return unit


@kamodofy
def rho(x=np.array([3, 4, 5])):
    """rho means density"""
    return x ** 2


@kamodofy
def divide(x=np.array([3, 4, 5])):
    """rho means density"""
    return x / 2


@np.vectorize
def myfunc(x, y, z):
    return x + y + z


def defaults_fun(var1, var2, default_var1=10, default_var2=20):
    return var1 + var2 + default_var1 + default_var2


@kamodofy
@gridify(x=np.linspace(-3, 3, 20),
         y=np.linspace(-5, 5, 10),
         z=np.linspace(-7, 7, 30))
def grid_fun(xvec):
    """my density function, returning array of size N

    xvec should have shape (N,3)"""
    return xvec[:, 0]


def test_kamodofy():
    comparison = rho.data == np.array([9, 16, 25])
    assert comparison.all()
    comparison = divide(np.array([6, 13, 2])) == np.array([3, 13 / 2, 1])
    assert comparison.all()


def test_sort_symbols():
    myexpr = parse_expr('x+rho+a+b+c')
    assert sort_symbols(myexpr.free_symbols) == (Symbol('a'), Symbol('b'), Symbol('c'), Symbol('rho'), Symbol('x'))


def test_valid_args():
    assert valid_args(myfunc, dict(x='a', other='you')) == OrderedDict([('x', 'a')])
    assert valid_args(myfunc, dict(y='b', other='you')) == OrderedDict([('y', 'b')])
    assert valid_args(myfunc, dict(z='c', other='you')) == OrderedDict([('z', 'c')])
    assert valid_args(myfunc, dict(x='a', y='b', z='c', other='you')) == OrderedDict(
        [('x', 'a'), ('y', 'b'), ('z', 'c')])


def test_eval_func():
    assert eval_func(myfunc, dict(x='a', y='b', z='c', other='you')) == myfunc('a', 'b', 'c')
    assert eval_func(myfunc, dict(x=1, y=2, z=3, other=100)) == myfunc(1, 2, 3)


def test_get_defaults():
    assert get_defaults(defaults_fun) == {'default_var1': 10, 'default_var2': 20}


def test_cast_0_dim():
    to = np.array([1, 2, 3])
    casted_array = cast_0_dim(np.array(1), to=to)
    assert casted_array.shape == to.shape


def test_concat_solution():
    solutions = concat_solution((kamodofy(lambda x=np.array([1, 2, 3]): x ** 2) for i in range(1, 4)), 'y')
    expected = [1, 4, 9, np.nan, 1, 4, 9, np.nan, 1, 4, 9, np.nan, ]
    assert ((solutions['y'] == expected) | (numpy.isnan(solutions['y']) & numpy.isnan(expected))).all()


def test_get_unit_quantity():
    from kamodo import get_unit
    mykm = get_unit_quantity('mykm', 'km', scale_factor=2)
    mygm = get_unit_quantity('mygm', 'gram', scale_factor=4)
    assert str(mykm.convert_to(get_unit('m'))) == '2000*meter'
    assert str(mygm.convert_to(get_unit('kg'))) == 'kilogram/250'


def test_substr_replace():
    mystr = "replace this string"
    mystr = substr_replace(mystr, [
        ('this', 'is'),
        ('replace', 'this'),
        ('string', 'replaced')
    ])
    assert mystr == 'this is replaced'


def test_beautify_latex():
    beautified = beautify_latex('LEFTx__plus1RIGHTminusLEFTacommabRIGHT')
    assert beautified == '\\left (x__+1\\right )-\\left (a,b\\right )'


def test_arg_to_latex():
    my_latex = arg_to_latex('x__iplus1')
    assert my_latex == 'x^{i+1}'


def test_valid_args():
    def f(x, y, z):
        return x + y + z

    args = valid_args(f, dict(x=1, y=2, z=3, g=4))
    assert args['x'] == 1
    assert ('g' not in args)


def test_eval_func():
    def f(x, y):
        return x + y

    assert eval_func(f, dict(x=1, y=2, z=3)) == 1 + 2


def test_simulate():
    def update_y(x):
        return x + 1

    def update_x(y):
        return y - 2

    state_funcs = OrderedDict([
        ('y', update_y),
        ('x', update_x),
        ('t', lambda t, dt: t + dt)
    ])

    simulation = simulate(state_funcs,
                          x=3,  # initial conditions
                          t=0,
                          dt=1,
                          steps=10)

    for state in simulation:
        pass

    assert state['x'] == -7
    assert state['y'] == -5
    assert state['t'] == 10


def test_meta_units():
    @kamodofy(units='kg')
    def mass(x):
        return x

    assert mass.meta['units'] == 'kg'


def test_meta_data():
    @kamodofy(units='kg')
    def mass(x=list(range(5)), y=3):
        return [x_ + y for x_ in x]

    assert mass.data[-1] == 7


def test_meta_data_kwargs():
    @kamodofy(x=3, y=3)
    def e(x, y):
        return x + y

    assert e.data == 6

    @kamodofy(x=3)
    def f(x=5, y=3):
        return x + y

    assert f.data == 6

    @kamodofy()
    def g(x, y):
        return x + y

    assert g.data == None
    assert g.meta['units'] == ''

    @kamodofy(data=3)
    def h(x, y):
        return x + y

    assert h.data == 3


def test_repr_latex():
    @kamodofy(equation='$f(x) = x_i^2$')
    def f(x):
        return x

    assert f._repr_latex_() == 'f(x) = x_i^2'

    @kamodofy
    def g(x, y, z):
        return x
    print(g._repr_latex_())

    assert g._repr_latex_() == r'$g{\left(x,y,z \right)} = \lambda{\left(x,y,z \right)}$'


def test_bibtex():
    bibtex = """
@phdthesis{phdthesis,
  author       = {Peter Joslin},
  title        = {The title of the work},
  school       = {The school of the thesis},
  year         = 1993,
  address      = {The address of the publisher},
  month        = 7,
  note         = {An optional note}
}"""

    @kamodofy(citation=bibtex)
    def h(x):
        '''This equation came out of my own brain'''
        return x

    assert '@phdthesis' in h.meta['citation']


def test_gridify():
    from kamodo import Kamodo
    kamodo = Kamodo(grids=grid_fun)
    assert kamodo.grids().shape == (10, 20, 30)


def test_symbolic():
    assert symbolic(['a']) == UndefinedFunction('a')
    assert symbolic(['a', 'b']) == [UndefinedFunction('a'), UndefinedFunction('b')]


def test_pad_nan():
    array = np.array([[1, 2], [1, 2]])
    solution = pad_nan(array)
    expected = np.array([[1, 2], [1, 2], [np.nan, np.nan]])
    comparison = ((solution == expected) | (numpy.isnan(solution) & numpy.isnan(expected)))
    assert comparison.all()

    array = np.array([1, 2])
    solution = pad_nan(array)
    expected = np.array([1, 2, np.nan])
    comparison = ((solution == expected) | (numpy.isnan(solution) & numpy.isnan(expected)))
    assert comparison.all()
