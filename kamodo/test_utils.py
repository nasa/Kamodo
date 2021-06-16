from collections import OrderedDict

import numpy as np
import pytest
from sympy import sympify as parse_expr, Symbol, symbols
from sympy.core.function import UndefinedFunction

from kamodo import Kamodo
from kamodo.util import kamodofy, gridify, sort_symbols, valid_args, eval_func, get_defaults, cast_0_dim, \
    concat_solution, get_unit_quantity, substr_replace, beautify_latex, arg_to_latex, simulate, pad_nan, \
    pointlike, solve, event, is_function

from .util import serialize, deserialize
import pandas as pd
import json
from .util import LambdaGenerator
from .util import curry, partial, get_args

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


@pointlike(squeeze='0')
def squeeze_error_point(x=np.linspace(0, 8 * np.pi, 100)):
    return np.sin(x)


@pointlike(squeeze=0)
def squeeze_point(x=np.linspace(0, 8 * np.pi, 100)):
    return np.sin(x)


@pointlike()
def points(x=np.linspace(0, 8 * np.pi, 100)):
    return np.sin(x)


def fprime(x):
    return x * x


@event
def boundry(x, s=np.array([1, 2])):
    r = np.linalg.norm(s)

    if np.isnan(r):
        result = 0
    else:
        result = r - 1
    return result


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
    def f(x):
        return x + 'a'

    with pytest.raises(TypeError) as error:
        eval_func(f, dict(x=1))

    assert eval_func(myfunc, dict(x='a', y='b', z='c', other='you')) == myfunc('a', 'b', 'c')
    assert eval_func(myfunc, dict(x=1, y=2, z=3, other=100)) == myfunc(1, 2, 3)
    assert 'unsupported operand type(s) for +' in str(error)


def test_get_defaults():
    assert get_defaults(defaults_fun) == {'default_var1': 10, 'default_var2': 20}


def test_cast_0_dim():
    to = np.array([1, 2, 3])
    casted_array = cast_0_dim(np.array(1), to=to)
    assert casted_array.shape == to.shape


def test_concat_solution():
    solutions = concat_solution((kamodofy(lambda x=np.array([1, 2, 3]): x ** 2) for i in range(1, 4)), 'y')
    expected = [1, 4, 9, np.nan, 1, 4, 9, np.nan, 1, 4, 9, np.nan, ]
    assert ((solutions['y'] == expected) | (np.isnan(solutions['y']) & np.isnan(expected))).all()


def test_get_unit_quantity():
    from kamodo import get_unit
    mykm = get_unit_quantity('mykm', 'km', scale_factor=2)
    mygm = get_unit_quantity('mygm', 'gram', scale_factor=4)
    assert str(mykm.convert_to(get_unit('m'))) == '2000.0*meter'
    assert str(mygm.convert_to(get_unit('kg'))) == '0.004*kilogram'


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


def test_simulate_error():
    def update_y(x):
        return x + '1'

    def update_x(y):
        return y - 2

    state_funcs = OrderedDict([
        ('y', update_y),
        ('x', update_x),
        ('t', lambda t, dt: t + dt)
    ])

    with pytest.raises(TypeError) as error:
        simulation = simulate(state_funcs,
                              x=3,
                              t=0,
                              dt=1,
                              steps=10)

        for state in simulation:
            pass

    assert "y:unsupported operand type(s) for +" in str(error)


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



def test_pad_nan():
    array = np.array([[1, 2], [1, 2]])
    solution = pad_nan(array)
    expected = np.array([[1, 2], [1, 2], [np.nan, np.nan]])
    comparison = ((solution == expected) | (np.isnan(solution) & np.isnan(expected)))
    assert comparison.all()

    array = np.array([1, 2])
    solution = pad_nan(array)
    expected = np.array([1, 2, np.nan])
    comparison = ((solution == expected) | (np.isnan(solution) & np.isnan(expected)))
    assert comparison.all()

    with pytest.raises(NotImplementedError) as error:
        array = np.array([[[1], [2], [3]], [[1], [2], [3]]])
        pad_nan(array)

def test_pointlike():
    assert Kamodo(points=squeeze_point).points().shape == (100,)
    assert Kamodo(points=points).points().shape == (1, 100)
    assert Kamodo(points=squeeze_error_point).points().shape == (1, 100)


def test_solve():
    kamodo = Kamodo()
    seeds = np.array([1, 2])
    kamodo['fprime'] = solve(fprime, seeds, 'x', (0, 30), events=boundry)

    assert kamodo.fprime().shape == (2, 2)


def test_is_function():
    assert is_function(parse_expr('f(g(x))'))
    assert is_function(symbols('g', cls=UndefinedFunction))
    assert not is_function(symbols('x'))

def test_serialize_np():
    x = np.linspace(-5, 5, 12).reshape((3, 4))

    assert deserialize(serialize(x)).shape == (3, 4)
    assert deserialize(serialize(x))[0, 0] == -5


    x_json = json.dumps(x, default=serialize)
    x_ = json.loads(x_json, object_hook=deserialize)
    assert (x_ == x).all()

def test_serialize_pd():
    t = pd.date_range('Jan 1, 2021', 'Jan 11, 2021', freq='H')
    
    t_json = json.dumps(t, default=serialize)
    t_ = json.loads(t_json, object_hook=deserialize)
    assert (t_ == t).all()


    s = pd.Series(np.linspace(-5, 5, 10), index=list(range(10)))
    s_json = json.dumps(s, default=serialize)
    s_ = json.loads(s_json, object_hook=deserialize)
    assert (s_ == s).all()

    df = pd.DataFrame(np.linspace(-5, 5, 10), index=list(range(10)))
    df_json = json.dumps(df, default=serialize)
    df_ = json.loads(df_json, object_hook=deserialize)
    assert (df_ == df).all().all()

def test_serialize_generator():
    gen = (lambda x=b: x**2 for b in np.linspace(-5,5,12).reshape((4,3)))
    # need to reproduce the generator for comparison
    gen_duplicate = (lambda x=b: x**2 for b in np.linspace(-5,5,12).reshape((4,3)))
    gen_json = json.dumps(gen, default=serialize)
    gen_ = json.loads(gen_json, object_hook=deserialize)

    for f, f_ in zip(gen_, gen_duplicate):
        assert(f() == f_()).all()
        f_defaults = get_defaults(f)
        f__defaults = get_defaults(f_)
        for k, v in f_defaults.items():
            assert (v == f__defaults[k]).all()

def test_serialize_generator_pd():
    t = pd.date_range('Jan 1, 2021', 'Jan 11, 2021', freq='H')

    gent = (lambda t=t: (t - t0).total_seconds() for t0 in t[::10])

    for f in deserialize(serialize(gent)):
        assert len(f()) == len(t)
    assert f()[-1] == 0

def test_multiply_lambdagen():
    lamb = LambdaGenerator((lambda x=np.arange(i): x**2 for i in range(5)))
    lamb2 = LambdaGenerator((lambda x=np.arange(i): x**2 for i in range(5)))
    for r in lamb*lamb2:
        r()

def test_divide_lambdagen():
    lamb = LambdaGenerator(((lambda x=np.linspace(1, 5, 10): x**2) for i in range(5)))
    lamb2 = LambdaGenerator(((lambda x=np.linspace(1, 5, 10): x**2) for i in range(5)))
    for r in lamb/lamb2:
        r()

def test_get_args():
    def f(x, y, z):
        return x+y+z

    assert get_args(f) == ('x', 'y', 'z')

def test_curry():
    @curry
    def arimean(*args):
        return sum(args) / len(args)

    assert arimean(-2, -1, 0, 1, 2)() == 0
    assert arimean(-2)(-1)(0)(1)(2)() == 0

def test_partial_decorator():
    @partial(z=1, verbose=True)
    def f(x, y=2, z=5):
        """simple function"""
        return x + y + z
    assert f(2,3) == 2+3+1
    

def test_partial_order():
    @partial(y=1, verbose=True)
    def f(x, y=2, z=5):
        """simple function"""
        return x + y + z
    assert f(3,4) == 3+1+4

def test_partial_no_kwargs():
    @partial(verbose=True)
    def f(x, y=2, z=5):
        """simple function"""
        return x + y + z
    assert f(3,4) == 3+4+5


def test_partial_inline():
    def f(x, y=2, z=5):
        """simple function"""
        return x + y + z
    g = partial(f, y=5, verbose=True)
    assert(g(3,4) == 3+5+4)

def test_partial_no_defaults():
    @partial(x=1, y=2, z=3, verbose=True)
    def f(x, y, z):
        return x + y + z
    assert f() == 1 + 2 + 3

def test_partial_missing_defaults():
    @partial(x=1, y=2, verbose=True)
    def f(x,y,z):
        return x+y+z
    assert f(3) == 1 + 2 + 3

def test_partial_required_args():
    """need to make z a required argument"""
    @partial(x=1, y=2, verbose=True)
    def f(x,y,z):
        return x+y+z
    try:
        f()
    except TypeError as m:
        assert 'missing' in str(m)

def test_partial_docs():
    @partial(y=2)
    def f(x, y=3, z=4):
        """my docs"""
        return x + y + z
    
    assert 'my docs' in f.__doc__
    assert 'f(x, z=4) = f(x, y=2, z)' in f.__doc__

    @partial(y=2, z=1)
    def f(x, y=3, z=4):
        """my docs"""
        return x + y + z

    assert 'Calling f(x, y, z) for fixed y, z' in f.__doc__


    