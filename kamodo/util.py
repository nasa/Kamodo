# -*- coding: utf-8 -*-
"""
Copyright © 2017 United States Government as represented by the Administrator, National Aeronautics and Space Administration.
No Copyright is claimed in the United States under Title 17, U.S. Code.  All Other Rights Reserved.
"""
import copy
import inspect
import json
import sys
import tempfile
import types
from collections import OrderedDict, defaultdict
from datetime import datetime
from inspect import getfullargspec

import forge
import numpy as np
import pandas as pd
import sympy
from decorator import decorate
from numpy.distutils.exec_command import exec_command
from scipy.integrate import solve_ivp
from sympy import Add, Mul, Pow, Tuple, sympify
from sympy import Function
from sympy import latex, Eq
from sympy import nsimplify
from sympy import symbols, Symbol
from sympy.core.compatibility import reduce, Iterable
from sympy.core.function import UndefinedFunction
from sympy.physics import units
from sympy.physics import units as sympy_units
from sympy.physics.units import Dimension
from sympy.physics.units import UnitSystem
from sympy.physics.units.quantities import Quantity
from sympy.physics.units.systems.si import dimsys_SI
from sympy.physics.units.util import _get_conversion_matrix_for_expr
from sympy.utilities.autowrap import ufuncify


def get_unit_quantity(name, base, scale_factor, abbrev=None, unit_system='SI'):
    '''Define a unit in terms of a base unit'''
    unit = units.Quantity(name, abbrev=abbrev)
    base_unit = getattr(sympy_units, base)

    try:
        # sympy >= 1.5 raises the following warning:
        #   Use unit_system.set_quantity_dimension or
        # <unit>.set_global_relative_scale_factor
        unit.set_global_relative_scale_factor(scale_factor, base_unit)
    except AttributeError:
        unit.set_dimension(base_unit.dimension)
        unit.set_scale_factor(scale_factor * base_unit, unit_system=unit_system)

    return unit


unit_subs = dict(
    nT=get_unit_quantity('nanotesla', 'tesla', .000000001, 'nT', 'SI'),
    R_E=get_unit_quantity('earth radii', 'm', 6.371e6, 'R_E', 'SI'),
    R_S=get_unit_quantity('solar radii', 'm', 6.957e8, 'R_S', 'SI'),
    erg=get_unit_quantity('erg', 'J', .0000001, 'erg', 'SI'),
    nPa=get_unit_quantity('nanopascals', 'pascal', .000000001, 'nPa', 'SI'),
    cc=sympy_units.cm ** 3,
    AU=get_unit_quantity('astronomical unit', 'm', 1.496e+11, 'AU', 'SI'),
    arcsec=get_unit_quantity('arc seconds', 'degrees', 1. / 3600, '\"', 'SI'),
    hr=get_unit_quantity('hour', 's', 3600., 'hr', 'SI')
    # TECU = get_unit_quantity('TECU','1/m**2',10**16, 'TECU','SI')
)

sympy_units.erg = unit_subs['erg']

prefix_dict = sympy_units.prefixes.PREFIXES  # built-in dictionary of Prefix instances
# test_unit_subs={}  #dictionary to replace subs in kamodo.get_unit_quantities()
unit_list = ['m', 's', 'g', 'A', 'K', 'radian', 'sr', 'cd', 'mole', 'eV', 'Pa',
             'F', 'N',
             'V', 'Hz', 'C', 'W', 'Wb', 'H', 'S', 'Bq', 'Gy', 'erg', 'T']

# list of SI units included in sympy (likely not complete)

lambda_ = symbols('lambda', cls=UndefinedFunction)

for item in unit_list:
    unit_item = getattr(sympy_units, item)
    unit_subs[item] = unit_item
    for key in prefix_dict.keys():
        unit_subs[key + item] = get_unit_quantity(str(prefix_dict[key].name) +
                                                  str(unit_item.name),
                                                  str(unit_item.name),
                                                  float(prefix_dict[
                                                            key].scale_factor),
                                                  abbrev=key + item)

reserved_names = dir(sympy)


def get_kamodo_unit_system():
    """Same as SI but supports anglular frequency"""

    radian = sympy.physics.units.radian
    degree = sympy.physics.units.degree
    si_unit_system = UnitSystem.get_unit_system('SI')
    si_dimension_system = si_unit_system.get_dimension_system()

    angle = Dimension('angle', 'A')

    kamodo_dims = si_dimension_system.extend(
        new_base_dims=(angle,),
        new_derived_dims=[Dimension('angular_velocity')],
        new_dim_deps={Symbol('angular_velocity'): {Symbol('angle'): 1,
                                                   Symbol('time'): -1}})

    kamodo_units = si_unit_system.extend(
        (radian,),
        (radian, degree),
        dimension_system=kamodo_dims)

    return kamodo_units


kamodo_unit_system = get_kamodo_unit_system()


def get_ufunc(expr, variable_map):
    """Numerically optimize expression"""

    expr = sympify(expr, locals=variable_map)
    func = ufuncify(expr.free_symbols, expr)
    formula = 'f{} = {}'.format(tuple(expr.free_symbols), expr)
    return func, formula


def compile_fortran(source, module_name, extra_args='', folder='./'):
    """compile fortran source code"""
    with tempfile.NamedTemporaryFile('w', suffix='.f90') as fortran_file:
        fortran_file.write(source)
        fortran_file.flush()

        args = ' -c -m {} {} {}'.format(module_name, fortran_file.name,
                                        extra_args)
        command = 'cd "{}" && "{}" -c "import numpy.f2py as f2py;f2py.main()" {}'.format(
            folder, sys.executable, args)
        status, output = exec_command(command)
        return status, output, command


def substr_replace(name, name_maps):
    """replaces all substrings in name with those given by name_maps"""
    for old, new in name_maps:
        name = name.replace(old, new)
    return name


# This should be in yaml
def beautify_latex(expr):
    """convert string to latex-compatible expression"""
    return substr_replace(expr, [
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


def kamodofy(
        _func=None,
        units='',
        arg_units=None,
        data=None,
        update=None,
        equation=None,
        citation=None,
        hidden_args=[],
        **kwargs):
    """
    Adds `meta` and `data` attributes to functions for registering with Komodo objects.
    
    ** inputs **:

    * _func: function to wrap
    * units: (optional) physical output units
    * arg_units: (optional) dictionary { arg : str unit} containing physical input units
    * data: if supplied, set f.data = data, if not supplied, set f.data = f(), assuming it can be called with no arguments.
      If f cannot be called with no arguments, will set f.data = None
    * update: name of another function's argument to update (see [simulation api](../notebooks/Kamodo/#simulation-api))
    * equation: str representing right-hand-side of the function
    * citation: str reference for publication
    * hidden_args: arguments of function to hide from latex rendering
    * kwargs: other key word arguments

    ** returns **: the decorated function with the following attributes

    * meta is a dictionary containing
        * units: physical output units (str)
        * arg_units: dictionary { arg : str unit}
        * equation: latex str representing right-hand-side of the function
        * citation: str reference for publication
        * hidden: str list of arguments to hide from latex rendering
    * data: default data representing expected function output for default arguments
    * update: name of another function's argument to update


    ** usage **:

    ```python
    @kamodofy(units='kg/cm^2', arg_units=dict(x='cm'), citation='Pembroke et. al 2022', hidden_args=['verbose'])
    def myfunc(x=30, verbose=True):
        return x**2
    myfunc.meta
    ```
    ```console
    {'units': 'kg/cm^2',
     'arg_units': {'x': 'cm'},
     'citation': 'Pembroke et. al 2022',
     'equation': None,
     'hidden_args': ['verbose']}
    ```
    The above metadata is used by Kamodo objects for function registration. Similarly, a `data` attribute is attached which represents the output of the function when called with no arguments:

    ```python
    myfunc.data
    ```
    ```console
    900
    ```

    """

    def decorator_kamodofy(f):
        f.meta = dict(
            units=units,
            arg_units=arg_units,
            citation=citation,
            equation=equation,
            hidden_args=hidden_args)

        if citation is not None:
            if f.__doc__ is not None:
                f.__doc__ = f.__doc__ + '\n\ncitation: {}'.format(citation)
            else:
                f.__doc__ = 'citation: {}'.format(citation)
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
            # lambda_ = symbols('lambda', cls=UndefinedFunction)
            lambda_ = Function('lambda')
            latex_eq = latex(Eq(lhs, lambda_(*lhs.args)))
            f._repr_latex_ = lambda: "${}$".format(latex(latex_eq))

        return decorate(f, decorator_wrapper)  # preserves signature

    if _func is None:
        return decorator_kamodofy
    else:
        return decorator_kamodofy(_func)


def sort_symbols(symbols):
    symbols_ = list(symbols)
    symbols_.sort(key=str)
    return tuple(symbols_)


def valid_args(f, kwargs):
    '''Extract arguments from kwargs that appear in f'''
    valid = OrderedDict()
    if type(f) == np.vectorize:
        for a in inspect.getfullargspec(f.pyfunc).args:
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


def get_args(func):
    sig = inspect.signature(func)
    return tuple(sig.parameters.keys())


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
    if isinstance(array, types.GeneratorType):
        return np.array(list(array))
    try:
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
    except:
        raise NotImplementedError(
            'cannot handle {} of type {}'.format(array, type(array)))


def concat_solution(gen, variable):
    """combine solutions of a function generator
    iterates through a function generator and extracts defaults for each function
    these solutions are then padded with nans and stacked
    stacking is vertically for N-d, horizontally for 1-d
    """
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
existing_plot_types.index.set_names(['nargs', 'arg shapes', 'out shape'],
                                    inplace=True)
existing_plot_types.columns = ['Plot Type', 'notes']


def gridify(_func=None, order='A', squeeze=True, **defaults):
    """Given a function _func(xvec) taking a single variable of shape (n,dim) and defaults
    (e.g. x(L), y(M), z(N)) `gridify` returns a new function `_func(x,y,z)` that calls `_func` with points 
    `n=LxMxN`, reshaping the result to `(M, L, N)` (`order='A'`, default) or `(L, M, N)` (`order='C'`).
    See [np.meshgrid](https://numpy.org/doc/stable/reference/generated/numpy.meshgrid.html)
    and [np.reshape](https://numpy.org/doc/stable/reference/generated/numpy.reshape.html).

    ** inputs **:

    * _func: kamodo function
    * order: 'A' str (default) passed to reshape.
        * order = 'A': use indexing='xy' in meshgrid
        * order = 'C': use indexing='ij' in meshgrid
    * squeeze: True (default) passed to reshape before returning

    ** returns **: mutated function

    ** usage **:

    Conceptually, `@gridify` converts point-like interpolators to grid-based interpolators.

    Suppose we have a function `r(xvec)` that takes an array of shape `(n,2)` and returns the magnitude `r` of each point.
    By applying `@gridify`, we convert `r(xvec)` to `r(x,y)`:

    ```python
    @gridify(x=np.linspace(-3,3,5), y=np.linspace(-5,5,11), order='A')
    def r(xvec):
        return np.linalg.norm(xvec, axis=-1)
    ```

    The result is automatically reshaped to match the input arrays for `x,y`:
    
    ```python
    r(x=[2,3], y=[3,4,5]).shape # (3, 2)
    ```
    
    In addition, `r` receives defaults for `x` and `y`:

    ```python
    r().shape # (11,5)
    ```

    To see the defaults, use the `get_defaults` function

    ```py

    from kamodo import get_defaults

    defaults = get_defaults(r) 
    defaults['x'] # (5,)
    defaults['y'] # (11,)
    ```
    
    We can generate "slices" for fixed values of `x` or `y`:

    ```python
    r(x=0) # array([5., 4., 3., 2., 1., 0., 1., 2., 3., 4., 5.])
    r(y=0) # array([3. , 1.5, 0. , 1.5, 3. ])
    ```

    Use `order` to control the shapes of returned arrays (row vs column major).

    ```python
    @gridify(x=np.linspace(-3,3,5), y=np.linspace(-5,5,11), z=np.linspace(-1,1,13), order='A')
    def r(xvec):
        return np.linalg.norm(xvec, axis=-1)

    r().shape # (11, 5, 13) corresponds to y, x, z

    @gridify(x=np.linspace(-3,3,5), y=np.linspace(-5,5,11), z=np.linspace(-1,1,13), order='C')
    def r(xvec):
        return np.linalg.norm(xvec, axis=-1)

    r().shape # (5, 11, 13) corresponds to x, y, z
    ```

    By default, the output array will be squeezed if a dimension has size 1.
    Use `squeeze=False` to disable this behavior:

    ```python
    @gridify(x=np.linspace(-3,3,5), y=np.linspace(-5,5,11), squeeze=False)
    def r(xvec):
        return np.linalg.norm(xvec, axis=-1)

    r(y=0).shape # (1, 5)
    r(x=0).shape # (11, 1)
    ```

    """

    def decorator_gridify(f):

        signature = []
        for arg, arg_default in defaults.items():
            signature.append(forge.arg(arg, default=arg_default))

        if order == 'A':
            indexing = 'xy'
        else:
            indexing = 'ij'

        @forge.sign(*signature)
        def wrapped(**kwargs):
            coordinates = np.meshgrid(
                *kwargs.values(),
                indexing=indexing,
                sparse=False,
                copy=False)
            points = np.column_stack([c.ravel() for c in coordinates])

            if squeeze:
                out_shape = [-1] + list(coordinates[0].shape)
                return np.squeeze(f(points).reshape(out_shape, order=order))
            else:
                out_shape = list(coordinates[0].shape)
                return f(points).reshape(out_shape, order=order)

        wrapped.__name__ = f.__name__
        wrapped.__doc__ = f.__doc__

        return decorate(wrapped, decorator_wrapper)

    if _func is None:
        return decorator_gridify
    else:
        return decorator_gridify(_func)


def pointlike(_func=None, signature=None, otypes=[float], squeeze=None):
    """Transforms a single-argument function to one that accepts m points of dimension n
    pointlike wraps [np.vectorize](https://numpy.org/doc/stable/reference/generated/numpy.vectorize.html)

    ** inputs **:

    * _func: kamodo function
    * signature: (optional) Generalized universal function signature, e.g., (m,n),(n)->(m) for vectorized matrix-vector multiplication. 
    * otypes: list(types) default is [float] The output data type. It must be specified as either a string of typecode characters or a list of data type specifiers. There should be one data type specifier for each output.
    * squeeze: axis on which to squeeze result 

    ** returns **: modified function 

    """

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
                    return np.vectorize(f, otypes=otypes, signature=signature)(
                        *args, **kwargs).squeeze(squeeze)
                except:
                    return np.vectorize(f, otypes=otypes, signature=signature)(
                        *args, **kwargs)
            else:
                return np.vectorize(f, otypes=otypes, signature=signature)(
                    *args, **kwargs)

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
                seed_numbers = np.ones(
                    len(interval_bounded)) * i  # *len(directions) + 1*(d > 0)
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

            isolution = np.floor(s.real).astype(int) * len(directions) + (
                    s.imag > 0)

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
            return pd.DataFrame(np.vstack(results),
                                index=index_).drop_duplicates()

        solution.__name__ = varname

        return decorate(solution, decorator_wrapper)  # preserves signature

    if fprime is None:
        return decorator_solve
    else:
        return decorator_solve(fprime)


def convert_unit_to(expr, target_units, unit_system=kamodo_unit_system,
                    raise_errors=True):
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
        return Eq(convert_unit_to(expr.lhs, target_units, unit_system),
                  convert_unit_to(expr.rhs, target_units, unit_system))
    # if type(type(expr)) is UndefinedFunction:
    if is_function(expr):
        # print('undefined input expr:{}'.format(expr))
        return nsimplify(expr, rational=True)

    if isinstance(expr, Add):
        return Add.fromiter(
            convert_unit_to(i, target_units, unit_system) for i in expr.args)

    expr = sympify(expr)

    if not isinstance(expr, Quantity) and expr.has(Quantity):
        try:
            expr = expr.replace(
                lambda x: isinstance(x, Quantity),
                lambda x: x.convert_to(target_units, unit_system))
        except OSError:
            raise OSError('problem converting {} to {}\n{}'.format(
                expr, target_units, unit_system))

    def get_total_scale_factor(expr):
        if isinstance(expr, Mul):
            return reduce(lambda x, y: x * y,
                          [get_total_scale_factor(i) for i in expr.args])
        elif isinstance(expr, Pow):
            return get_total_scale_factor(expr.base) ** expr.exp
        elif isinstance(expr, Quantity):
            return unit_system.get_quantity_scale_factor(expr)
        return expr

    if expr == target_units:
        return expr

    depmat = _get_conversion_matrix_for_expr(expr, target_units, unit_system)
    if depmat is None:
        if raise_errors:
            raise NameError(
                'cannot convert {} to {} {}'.format(expr, target_units,
                                                    unit_system))
        return nsimplify(expr, rational=True)

    expr_scale_factor = get_total_scale_factor(expr)
    result = expr_scale_factor * Mul.fromiter(
        (1 / get_total_scale_factor(u) * u) ** p for u, p in
        zip(target_units, depmat))
    return nsimplify(result, rational=True)


def get_expr_unit(expr, unit_registry, verbose=False):
    '''Get units from an expression'''
    if is_function(expr):
        for func in unit_registry:
            if type(expr) == type(func):
                # b(a) = b(x)
                if verbose:
                    print('get_expr_unit: found match {} for {}'.format(func,
                                                                        expr))
                # {f(x): f(cm), f(cm): kg}
                func_units = unit_registry[func]
                if func_units in unit_registry:
                    return unit_registry[func_units]
                return func_units
        if verbose:
            print('get_expr_unit: no match found for {}'.format(expr))

    if isinstance(expr, Mul):
        results = []
        for arg in expr.args:
            result = get_expr_unit(arg, unit_registry, verbose)
            if result is not None:
                results.append(result)
        if len(results) > 0:
            return Mul.fromiter(results)
        else:
            return None

    if len(unit_registry) > 0:
        expr_unit = expr.subs(
            unit_registry, simultaneous=False).subs(
            unit_registry, simultaneous=False)
    else:
        expr_unit = expr

    if len(expr_unit.free_symbols) > 0:
        return None

    if isinstance(expr_unit, Add):
        # use the first term
        arg_0 = expr_unit.args[0]
        convert_unit_to(expr_unit, arg_0, kamodo_unit_system)
        result = arg_0
    else:
        result = expr_unit

    # if len(result.free_symbols) > 0:
    #     return None

    return get_base_unit(result)


def get_arg_units(expr, unit_registry, arg_units=None):
    """retrieves units of expression arguments
    """
    if arg_units is None:
        arg_units = dict()

    if expr in unit_registry:
        # unit_registry: {f(cm,km): kg}
        if is_function(unit_registry[expr]):
            for arg_, arg_unit in zip(expr.args, unit_registry[expr].args):
                arg_units[arg_] = arg_unit
        return arg_units

    # 2*a(x)
    for arg in expr.args:
        arg_units = get_arg_units(arg, unit_registry, arg_units)
    return arg_units


def replace_args(expr, from_map, to_map):
    # func_symbol = type(expr)
    arg_map = dict()
    for arg in expr.free_symbols:
        if (arg in from_map) & (arg in to_map):
            from_unit = from_map[arg]
            to_unit = to_map[arg]
            try:
                assert get_dimensions(from_unit) == get_dimensions(to_unit)
                arg_map[arg] = convert_unit_to(
                    arg * to_unit,
                    from_unit,
                    kamodo_unit_system) / from_unit
            except:
                raise NameError(
                    'cannot convert from {} to {}'.format(from_unit, to_unit))
    return expr.subs(arg_map)


def unify_args(expr, unit_registry, to_symbol, verbose):
    """replaces arguments in an expression"""
    expr_units = get_arg_units(expr, unit_registry)
    to_units = get_arg_units(to_symbol, unit_registry)
    if verbose:
        print('unify_args: {} arg units {}'.format(expr, expr_units))
        print('unify_args: to arg units', to_units)
    expr = replace_args(expr, expr_units, to_units)
    if verbose:
        print('unify_args: converted expression:', expr)
    return expr


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
    if isinstance(expr, Add):
        if verbose:
            print('unify: Adding expression: {} -> {}'.format(
                expr, get_expr_unit(to_symbol, unit_registry, verbose)))
        return Add.fromiter(
            [unify(arg, unit_registry, to_symbol, verbose) for arg in
             expr.args])

    expr_unit = get_expr_unit(expr, unit_registry, verbose)

    if verbose:
        print('unify: {} unit {}'.format(expr, expr_unit))
        print('unify: {} symbols {}'.format(expr, expr.free_symbols))
        print('unify: {} symbols {}'.format(to_symbol, to_symbol.free_symbols))
    try:
        assert expr.free_symbols.issubset(to_symbol.free_symbols)
    except:
        raise NameError("{} arguments not in {}".format(
            expr.free_symbols, to_symbol.free_symbols))

    if is_function(expr):
        if verbose:
            print('unify: function expression: {}'.format(expr))
        if to_symbol is not None:
            if verbose:
                print('\nunify: to_symbol args: {}'.format(to_symbol.args))
                print('unify: to_symbol free symbols: {}'.format(
                    to_symbol.free_symbols))
                print('unify: expr args: {}'.format(expr.args))
                print('unify: expr free symbols: {}'.format(expr.free_symbols))

            for k, v in unit_registry.items():
                if isinstance(expr, type(k)):
                    if len(k.free_symbols) > 0:
                        if verbose:
                            print('unify: found matching {} -> {}'.format(expr,
                                                                          k))
                        arg_units = get_arg_units(k, unit_registry)
                        if verbose:
                            print('unify: func units:', arg_units)
                            print('unify: {}->{}'.format(expr.args, k.args))
                        expr_units = {}
                        for arg, sym in zip(expr.args, k.args):
                            to_unit = arg_units.get(sym)
                            from_unit = get_expr_unit(arg, unit_registry)
                            if (from_unit is not None) and (
                                    to_unit is not None):
                                expr_units[arg] = convert_unit_to(
                                    arg * from_unit,
                                    to_unit,
                                    kamodo_unit_system) / to_unit
                        expr = expr.subs(expr_units)

                        if verbose:
                            print('unify: replaced args', expr)

    expr = unify_args(expr, unit_registry, to_symbol, verbose)

    if (to_symbol is not None) & (expr_unit is not None):
        to_unit = get_expr_unit(to_symbol, unit_registry, verbose)
        if verbose:
            print('unify: to_unit {}'.format(to_unit))
        expr_dimensions = get_dimensions(expr_unit)
        to_dimensions = get_dimensions(to_unit)
        if expr_dimensions.compare(to_dimensions) == 0:
            if verbose:
                print('unify: {} [{}] -> to_symbol: {}[{}]'.format(
                    expr, expr_unit, to_symbol, to_unit))
            expr = convert_unit_to(
                expr * expr_unit,
                to_unit,
                kamodo_unit_system) / to_unit
        else:
            if verbose:
                print('unify: registry:')
                for k, v in unit_registry.items():
                    print('unify:\t{} -> {}'.format(k, v))
                print(
                    'compare:{}'.format(expr_dimensions.compare(to_dimensions)))

            error_msg = 'cannot convert {} [{}] {} to {}[{}] {}'.format(
                expr, expr_unit, expr_dimensions,
                to_symbol, to_unit, to_dimensions)
            raise NameError(error_msg)
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


def base_dimensions(d):
    dependencies = dimsys_SI.get_dimensional_dependencies(d)
    return Mul.fromiter(
        Pow(Dimension(base), exp_) for base, exp_ in dependencies.items())


def get_dimensions(unit):
    """get the set of basis units"""
    if hasattr(unit, 'dimension'):
        base_dims = base_dimensions(unit.dimension)
        if len(base_dims.args) == 1:
            return base_dims.args[0]
        else:
            return Mul.fromiter(
                [arg.args[0] for arg in base_dimensions(unit.dimension).args])
    if isinstance(unit, Mul):
        terms = [get_dimensions(arg) for arg in unit.args]
        return Mul.fromiter(terms)
    if isinstance(unit, Pow):
        base, exp = unit.as_base_exp()
        return Pow(get_dimensions(base), exp)
    return 1


def get_base_unit(expr):
    if hasattr(expr, 'dimension'):
        return expr
    if isinstance(expr, Mul):
        results = []
        for arg in expr.args:
            result = get_base_unit(arg)
            if result is not None:
                results.append(result)
        if len(results) > 0:
            return Mul.fromiter(results)
        else:
            return None
    if isinstance(expr, Pow):
        base, exp = expr.as_base_exp()
        return Pow(get_base_unit(base), exp)
    if isinstance(expr, Add):
        results = [get_dimensions(arg) for arg in expr.args]
        arg0_dimensions = results[0]
        for r in results:
            try:
                assert arg0_dimensions == r
            except:
                print(results)
                raise
        return get_base_unit(expr.args[0])
    return None


def is_function(expr):
    """returns True if expr is a function
    examples:
        f(x): True
        symbols('f', cls=UndefinedFunction): True
        x: False
    """
    if isinstance(expr, UndefinedFunction):
        return True
    return isinstance(type(expr), UndefinedFunction)


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


class NumpyArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NumpyArrayEncoder, self).default(obj)


def full_classname(o):
    module = o.__class__.__module__
    if module is None or module == str.__class__.__module__:
        return o.__class__.__name__
    return module + '.' + o.__class__.__name__


def serialize(obj):
    if isinstance(obj, (np.ndarray, np.generic)):
        if isinstance(obj, np.ndarray):
            return {
                # '__ndarray__': base64.b64encode(obj.tobytes()),
                '__ndarray__': obj.tolist(),
                'dtype': obj.dtype.str,
                'shape': obj.shape,
            }
        elif isinstance(obj, (np.bool_, np.number)):
            return {
                # '__npgeneric__': base64.b64encode(obj.tobytes()),
                '__ndarray__': obj.tolist(),
                'dtype': obj.dtype.str,
            }
    if isinstance(obj, pd.DataFrame):
        return {
            '__pddataframe__': obj.values.tolist(),
            '__index__': serialize(obj.index),
            'dtype': full_classname(obj),
        }
    if isinstance(obj, pd.Index):
        if isinstance(obj, pd.DatetimeIndex):
            return {
                '__datetime__': [_ for _ in map(datetime.isoformat, obj)],
                'dtype': full_classname(obj),
            }
        else:
            return {
                '__index__': obj.tolist(),
                'dtype': full_classname(obj),
            }
    if isinstance(obj, pd.Series):
        return {
            '__pdseries__': obj.values.tolist(),
            '__index__': serialize(obj.index),
            'dtype': full_classname(obj),
        }
    if isinstance(obj, set):
        return {'__set__': list(obj)}
    if isinstance(obj, tuple):
        return {'__tuple__': list(obj)}
    if isinstance(obj, complex):
        return {'__complex__': obj.__repr__()}

    if isinstance(obj, types.GeneratorType):
        return {'__lambdagen__': [{
            'params': {k: serialize(v) for k, v in get_defaults(func).items()},
            'result': serialize(func())} for func in obj]}

    if isinstance(obj, int):
        return obj

    # Let the base class default method raise the TypeError
    raise TypeError('Unable to serialise object of type {}'.format(type(obj)))


def lambdagen(obj):
    """create a generator of lambda functions"""
    for func_ in obj['__lambdagen__']:
        signature = []
        for arg, arg_default in func_['params'].items():
            signature.append(forge.arg(
                arg,
                default=deserialize(arg_default)))

        @forge.sign(*signature)
        def func(*args, **kwargs):
            """API function"""
            return deserialize(func_['result'])

        yield kamodofy(func)


def deserialize(obj):
    # convert obj into numpy, pandas
    if isinstance(obj, dict):
        if '__ndarray__' in obj:
            # return np.frombuffer(
            #     base64.b64decode(obj['__ndarray__']),
            #     dtype=np.dtype(obj['dtype'])
            # ).reshape(obj['shape'])
            return np.array(obj['__ndarray__'])
        if '__npgeneric__' in obj:
            # return np.frombuffer(
            #     base64.b64decode(obj['__npgeneric__']),
            #     dtype=np.dtype(obj['dtype'])
            # )[0]
            return np.array(obj['__npgeneric__'])
        if '__pddataframe__' in obj:
            return pd.DataFrame(
                obj['__pddataframe__'],
                index=deserialize(obj['__index__']))
        if '__pdseries__' in obj:
            return pd.Series(
                obj['__pdseries__'],
                index=deserialize(obj['__index__']))
        if '__datetime__' in obj:
            return pd.to_datetime(obj['__datetime__'])
        if '__index__' in obj:
            return pd.Index(obj['__index__'])
        if '__set__' in obj:
            return set(obj['__set__'])
        if '__tuple__' in obj:
            return tuple(obj['__tuple__'])
        if '__complex__' in obj:
            return complex(obj['__complex__'])
        if '__lambdagen__' in obj:
            return lambdagen(obj)

    return obj


# over-write the load(s)/dump(s) functions
def load(*args, **kwargs):
    kwargs['object_hook'] = deserialize
    return json.load(*args, **kwargs)


def loads(*args, **kwargs):
    kwargs['object_hook'] = deserialize
    return json.loads(*args, **kwargs)


def dump(*args, **kwargs):
    kwargs['default'] = serialize
    return json.dump(*args, **kwargs)


def dumps(*args, **kwargs):
    kwargs['default'] = serialize
    return json.dumps(*args, **kwargs)


def get_undefined_funcs(expr):
    """retrieve an expression's undefined functions"""
    return expr.atoms(sympy.function.AppliedUndef)


def sign_defaults(symbol, expr, composition):
    '''gets defaults from an expression using composition'''
    defaults = {}
    for f_ in get_undefined_funcs(expr):
        # includes {h(f(g)), f(g)}
        if str(f_) in composition:
            # ignores h(f(g))
            f_defaults = get_defaults(composition[str(f_)])
            # flatten defaults
            for arg, arg_default in f_defaults.items():
                defaults[arg] = arg_default
    arg_signatures = []
    # defaults have to go last, which may conflict with user's ordering
    symbol_args = list(symbol.args)
    for arg in symbol.args:
        str_arg = str(arg)
        if not (str_arg in defaults):
            arg_signatures.append(forge.arg(str_arg))
            symbol_args.remove(arg)
        else:
            continue

    for default_arg in symbol_args:
        str_arg = str(default_arg)
        arg_default = defaults[str(default_arg)]
        arg_signatures.append(forge.arg(str_arg, default=arg_default))
    # will raise an error if defaults are not last
    signature = forge.sign(*arg_signatures)
    return signature, defaults


class LambdaGenerator(object):
    def __init__(self, lambda_generator):
        """supports simple expressions for manipulating lambda generators"""
        self._lambda_generator = lambda_generator

    def __add__(self, other):
        if isinstance(other, LambdaGenerator):
            for func, gunc in zip(self._lambda_generator,
                                  other._lambda_generator):
                yield lambda: func() + gunc()
        else:
            for func in self._lambda_generator:
                yield lambda: func() + other

    def __sub__(self, other):
        if isinstance(other, LambdaGenerator):
            for func, gunc in zip(self._lambda_generator,
                                  other._lambda_generator):
                yield lambda: func() - gunc()
        else:
            for func in self._lambda_generator:
                yield lambda: func() - other

    def __mul__(self, other):
        if isinstance(other, LambdaGenerator):
            for func, gunc in zip(self._lambda_generator,
                                  other._lambda_generator):
                # what do we do with defaults?
                yield lambda: func() * gunc()
        else:
            for func in self._lambda_generator:
                yield lambda: func() * other

    def __truediv__(self, other):
        if isinstance(other, LambdaGenerator):
            for func, gunc in zip(self._lambda_generator,
                                  other._lambda_generator):
                # what do we do with defaults?
                yield lambda: func().__truediv__(gunc())
        else:
            for func in self._lambda_generator:
                yield lambda: func().__truediv__(other)

    def __floordiv__(self, other):
        if isinstance(other, LambdaGenerator):
            for func, gunc in zip(self._lambda_generator,
                                  other._lambda_generator):
                # what do we do with defaults?
                yield lambda: func().__floordiv__(gunc())
        else:
            for func in self._lambda_generator:
                yield lambda: func().__floordiv__(other)

    def __pow__(self, other):
        if isinstance(other, LambdaGenerator):
            for func, gunc in zip(self._lambda_generator,
                                  other._lambda_generator):
                yield lambda: func().__pow__(gunc())
        else:
            for func in self._lambda_generator:
                yield lambda: func().__pow__(other)


def get_bbox(fig):
    t0 = fig['data'][0]
    xmin = np.nanmin(t0['x'])
    xmax = np.nanmax(t0['x'])
    ymin = np.nanmin(t0['y'])
    ymax = np.nanmax(t0['y'])
    if 'z' in t0:
        zmin = np.nanmin(t0['z'])
        zmax = np.nanmax(t0['z'])
        for t in fig['data']:
            xmin = min(xmin, np.nanmin(t['x']))
            xmax = max(xmax, np.nanmax(t['x']))
            ymin = min(ymin, np.nanmin(t['y']))
            ymax = max(ymax, np.nanmax(t['y']))
            zmin = min(zmin, np.nanmin(t['z']))
            zmax = max(zmax, np.nanmax(t['z']))
        return xmin, xmax, ymin, ymax, zmin, zmax
    else:
        for t in fig['data']:
            xmin = min(xmin, np.nanmin(t['x']))
            xmax = max(xmax, np.nanmax(t['x']))
            ymin = min(ymin, np.nanmin(t['y']))
            ymax = max(ymax, np.nanmax(t['y']))
        return xmin, xmax, ymin, ymax


def set_aspect(fig):
    """sets aspect ratio of the scene based on bounding box of traces"""
    fig.layout.scene.aspectmode = 'manual'
    xmin, xmax, ymin, ymax, zmin, zmax = get_bbox(fig)
    fig.layout.scene.aspectratio = dict(x=xmax - xmin, y=ymax - ymin,
                                        z=zmax - zmin)
    return fig


def curry(func):
    """currying function
    borrowed from https://www.python-course.eu/currying_in_python.php
    """
    # to keep the name of the curried function:

    curry.__curried_func_name__ = func.__name__
    f_args, f_kwargs = [], {}

    def f(*args, **kwargs):
        """a curried function"""

        nonlocal f_args, f_kwargs
        if args or kwargs:
            f_args += args
            f_kwargs.update(kwargs)
            return f
        else:
            result = func(*f_args, **f_kwargs)
            f_args, f_kwargs = [], {}
            return result

    return f


def construct_signature(*args, **kwargs):
    """construct a signature
    usage:
        
        @forge.sign(*construct_signature('x','y',z=3))
        def f(*args, **kargs):
            pass
    """
    signature = []
    for arg in args:
        signature.append(forge.arg(arg))
    for k, v in kwargs.items():
        signature.append(forge.arg(k, default=v))
    return signature


def array_to_latex(arr, max_chars=3):
    """a latex-friendly render for arrays"""
    try:
        iter(arr)
    except:
        return arr

    if hasattr(arr, 'shape'):
        if len(arr.shape) > 1:
            return '[{},..]'.format(''.join(array_to_latex(arr[0], max_chars)))
        if len(arr) > max_chars:
            return '[{},..]'.format(','.join(arr[:max_chars].astype(str)))
        return '[{}]'.format(','.join(arr.astype(str)))

    if len(arr) > max_chars:
        return '[{},..]'.format(','.join([str(_) for _ in arr[:max_chars]]))
    else:
        return arr


def latex_repr_values(values_dict):
    sigs = []
    for k, v in values_dict.items():
        sigs.append('{} = {}'.format(k, getattr(v, '_repr_latex_',
                                                array_to_latex(v, 2))))
    return ','.join(sigs)


def partial(_func=None, **partial_kwargs):
    """A partial function decorator, reducing function signature to reflect partially assigned kwargs.
    
    ** inputs **:

    * _func: kamodo function

    * partial_kwargs: (dict) _func arguments to set

    ** returns **: updated function with reduced arguments

    ** usage **:

    ```python
    @partial(z=1)
    def f(x, y=2, z=5):
        return x + y + z
    assert f(2,3) == 2+3+1
    try:
        f(3,4,5)
    except TypeError as m:
        print(m) # wrapped() takes from 1 to 2 positional arguments but 3 were given
    ```

    Note: This decorator differs significantly from functools.partial in the following ways:

    * functools.partial updates the function defaults without actually eliminating arguments.
    * functools.partial raises TypeError when used as a @decorator
    """
    verbose = partial_kwargs.pop('verbose', False)

    def decorator_partial(f):
        orig_args = get_args(f)
        orig_defaults = get_defaults(f)
        orig_defaults.update(partial_kwargs)
        orig_latex_func = getattr(f,
                                  '_repr_latex_',
                                  lambda: '\\lambda ({})'.format(
                                      ','.join(orig_args)))
        orig_meta = copy.deepcopy(getattr(f, 'meta', {}))
        orig_equation = orig_meta.get('equation')
        if orig_equation is None:
            orig_equation = orig_latex_func()
        orig_meta['equation'] = orig_equation.strip("$")

        # remove partials from arg units dictionary
        orig_arg_units = orig_meta.get('arg_units')
        if orig_arg_units is not None:
            for _ in partial_kwargs:
                orig_arg_units.pop(_)
        orig_meta['arg_units'] = orig_arg_units

        if len(orig_defaults) > 0:
            orig_meta['equation'] += ', ' + latex_repr_values(partial_kwargs)
        new_latex_func = lambda: '${}$'.format(orig_meta['equation'])
        if verbose:
            print('partial kwargs', partial_kwargs)
            print('original args:', orig_args)
            print('new defaults', orig_defaults)
            print('original latex function', orig_latex_func)
            print('new latex function', new_latex_func)

        # collect only the arguments not assigned by partial
        sig_defaults = {}
        sig_args = []
        for arg in orig_args:
            if arg in partial_kwargs:
                continue
            if arg in orig_defaults:
                sig_defaults[arg] = orig_defaults[arg]
            else:
                sig_args.append(arg)

        if verbose:
            print('updated signature:', sig_args, sig_defaults)

        @forge.sign(*construct_signature(*sig_args, **sig_defaults))
        def wrapped(*args, **kwargs):
            """simple wrapper"""
            kwargs.update(partial_kwargs)
            if verbose:
                print('kwargs to pass:', kwargs)
            return f(*args, **kwargs)

        if verbose:
            print('wrapped docs', wrapped.__doc__)
        wrapped = decorate(wrapped, decorator_wrapper)

        wrapped.__name__ = f.__name__
        orig_docs = f.__doc__

        for k, v in sig_defaults.items():
            sig_args.append('{}={}'.format(k, v))
        doc_args = ', '.join(sig_args)
        rhs_args = []
        for arg in orig_args:
            if arg in partial_kwargs:
                rhs_args.append('{}={}'.format(arg, partial_kwargs[arg]))
            else:
                rhs_args.append(arg)

        doc_orig_args = ', '.join(rhs_args)
        wrapped.__doc__ = """Calling {fname}({orig_args}) for fixed {partial_keys}:
        {fname}({doc_args}) = {fname}({orig_args_fixed})\n""".format(
            orig_args_fixed=doc_orig_args,
            doc_args=doc_args,
            fname=f.__name__,
            orig_args=', '.join(orig_args),
            partial_keys=', '.join(partial_kwargs.keys()),
        )
        wrapped._repr_latex_ = new_latex_func
        wrapped.meta = orig_meta
        if orig_docs is not None:
            wrapped.__doc__ += "\n" + f.__doc__

        return wrapped

    if _func is None:
        return decorator_partial
    else:
        return decorator_partial(_func)
