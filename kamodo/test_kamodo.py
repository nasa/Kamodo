"""
Tests for kamodo.py

"""

import numpy as np
from sympy import symbols, Symbol
import pytest
from sympy.core.function import UndefinedFunction
import pandas as pd
from kamodo import Kamodo, get_unit, kamodofy, validate_units, Eq
import functools
from sympy import lambdify, sympify
from kamodo import get_abbrev
from .util import get_arg_units

def test_Kamodo_expr():
    a, b, c, x, y, z = symbols('a b c x y z')
    kamodo = Kamodo(a=x ** 2, verbose=True)
    try:
        assert kamodo.a(3) == 3 ** 2
    except:
        print(kamodo.symbol_registry)
        raise


def test_Kamodo_ambiguous_expr():
    a, b, c, x, y, z = symbols('a b c x y z')
    with pytest.raises(NameError):
        kamodo = Kamodo()
        kamodo['a(b,c)'] = x ** 2 + y ** 2
    kamodo = Kamodo()
    kamodo['a(x, y)'] = x ** 2 + y ** 2
    assert kamodo.a(3, 4) == 3 ** 2 + 4 ** 2


def test_Kamodo_callable():
    kamodo = Kamodo(f=lambda x: x ** 2, verbose=False)
    assert kamodo.f(3) == 9
    kamodo = Kamodo(f=lambda x, y: x ** 2 + y ** 3)
    assert kamodo.f(3, 4) == 3 ** 2 + 4 ** 3


def test_Kamodo_callable_array():
    kamodo = Kamodo(f=lambda x: x ** 2)
    assert (kamodo.f(np.linspace(0, 1, 107)) == np.linspace(0, 1, 107) ** 2).all()


def test_Kamodo_str():
    kamodo = Kamodo('f(a,x,b) = a*x+b')
    try:
        assert kamodo.f(3, 4, 5) == 3 * 4 + 5
    except:
        print(kamodo.signatures)
        print(list(kamodo.keys()))
        raise


def test_Kamodo_latex():
    kamodo = Kamodo('$f(a,x,b) = a^x+b $')
    assert kamodo.f(3, 4, 5) == 3 ** 4 + 5


def test_Kamodo_get():
    kamodo = Kamodo('$f(a,x,b) = a^x+b $')
    try:
        assert kamodo['f'](3, 4, 5) == 3 ** 4 + 5
    except:
        print(kamodo.signatures)
        print(list(kamodo.keys()))
        raise


def test_Kamodo_set():
    kamodo = Kamodo()
    kamodo['f(a,x,b)'] = '$a^x+b$'
    assert kamodo.f(2, 3, 4) == 2 ** 3 + 4


def test_Kamodo_mismatched_symbols():
    with pytest.raises(NameError):
        kamodo = Kamodo('$f(a,b) = a + b + c$', verbose=False)
        assert 'f' not in kamodo


def test_Kamodo_assign_by_expression():
    kamodo = Kamodo()
    f, x, y = symbols('f x y')
    kamodo['f(x, y)'] = x ** 2 + y ** 2
    assert kamodo['f'](3, 4) == 3 ** 2 + 4 ** 2
    assert kamodo.f(3, 4) == 3 ** 2 + 4 ** 2


def test_Kamodo_assign_by_callable():
    kamodo = Kamodo(f=lambda x: x ** 2, verbose=False)
    assert kamodo['f'](3) == 3 ** 2
    assert kamodo.f(3) == 3 ** 2
    Kamodo(f=lambda x, y: x + y)

    def rho(x, y):
        return x + y

    kamodo['g'] = rho
    assert kamodo.g(3, 4) == 3 + 4


def test_Kamodo_composition():
    # f, g = symbols('f g', cls=UndefinedFunction)
    # x = Symbol('x')
    kamodo = Kamodo(f="$x^2$", verbose=True)
    kamodo['g(x)'] = "$f(x) + x^2$"
    kamodo['h(x)'] = 'g**2 + x**2'
    
    try:
        assert kamodo.f(3) == 3 ** 2
        assert kamodo.g(3) == 3 ** 2 + 3 ** 2
        assert kamodo.h(3) == (3 ** 2 + 3 ** 2) ** 2 + 3 ** 2
    except:
        print(kamodo)
        raise


def test_Kamodo_reassignment():
    a, b, c, x, y, z, r = symbols('a b c x y z r')
    kamodo = Kamodo('f(x) = 3*x+5', verbose=True)
    kamodo['r(x,y)'] = x + y
    kamodo['r(x,y)'] = "3*x + y"
    assert kamodo.r(1, 1) != 2
    assert kamodo.r(1, 1) == 4


def test_multivariate_composition():
    kamodo = Kamodo(f='x**2', g=lambda y: y ** 3, verbose=True)
    kamodo['h(x,y)'] = 'f(x) + g(y)'
    kamodo['i(x,y)'] = 'x + f(h(x,y))'
    assert kamodo.h(3, 4) == 3 ** 2 + 4 ** 3
    assert kamodo.i(3, 4) == 3 + (3 ** 2 + 4 ** 3) ** 2


def test_special_numbers():
    kamodo = Kamodo()
    kamodo['f'] = "${}^x$".format(np.e)
    assert np.isclose(kamodo.f(1), np.e)


def test_function_registry():
    kamodo = Kamodo("f(x) = x**2", "g(y) = y**3")
    kamodo['r(x,y)'] = "(x**2 + y**2)**(1/2)"
    kamodo['h(x,y)'] = 'f + g + r'
    assert 'h(x, y)' in kamodo.signatures


def test_symbol_key():
    f = symbols('f', cls=UndefinedFunction)
    x = Symbol('x')
    f_ = f(x)
    kamodo = Kamodo(verbose=True)
    kamodo[f_] = x ** 2
    assert kamodo.f(3) == 9


def test_compact_variable_syntax():
    f = symbols('f', cls=UndefinedFunction)
    x = symbols('x')
    kamodo = Kamodo(f='x**2')  # f also registers f(x)
    kamodo['g(x)'] = 'f + x'
    assert kamodo.g(3) == 3 ** 2 + 3


def test_unit_registry():
    def rho(x, y):
        return x * y

    def v(x, y):
        return x + y

    v.meta = dict(units='km/s')
    kamodo = Kamodo(v=v, verbose=False)
    kamodo['rho[kg/cc]'] = rho
    try:
        assert kamodo.rho.meta['units'] == 'kg/cc'
        assert kamodo.v.meta['units'] == 'km/s'
    except:
        print('\n', pd.DataFrame(kamodo.signatures))
        raise

    with pytest.raises(NameError):
        kamodo['p[kg]'] = v
    assert 'p' not in kamodo


def test_to_latex():
    kamodo = Kamodo(f='x**2', verbose=True)
    assert str(kamodo.to_latex()) == r'\begin{equation}f{\left(x \right)} = x^{2}\end{equation}'
    kamodo = Kamodo(g='x', verbose=True)
    assert str(kamodo.to_latex()) == r'\begin{equation}g{\left(x \right)} = x\end{equation}'


def test_expr_conversion():
    kamodo = Kamodo('$a[km] = x$', verbose=True)
    print(kamodo.items())
    kamodo.a


def test_validate_units():
    f, x = symbols('f x')
    with pytest.raises(ValueError):
        validate_units(Eq(f, x), dict(f=get_unit('kg'), x=get_unit('m')))
    validate_units(Eq(f, x), dict(f=get_unit('kg'), x=get_unit('g')))

    lhs_units = validate_units(sympify('f(x)'), dict(f=get_unit('kg'), x=get_unit('m')))
    print(lhs_units)


def test_unit_conversion():
    kamodo = Kamodo('$a(x[m])[km/s] = x$',
                    '$b(y[cm])[m/s] = y$', verbose=True)
    kamodo['c(x[m],y[m])[m/s]'] = '$a + b$'
    assert kamodo.c(1, 2) == 1000 + 200


def test_get_unit():
    assert get_unit('kg/m^3') == get_unit('kg/m**3')


def test_unit_conversion_syntax():
    kamodo = Kamodo('rho[kg/m^3] = x', verbose=True)
    with pytest.raises(NameError):
        kamodo['d[kg]'] = 'rho'
        print(kamodo.detail())


def test_unit_composition():
    kamodo = Kamodo('m[kg] = x', verbose=True)
    kamodo['v[km/s]'] = 'y'
    kamodo['p(x,y)'] = 'm*v'
    try:
        assert get_unit(kamodo.signatures['p(x, y)']['units']) == get_unit('kg*km/s')
    except:
        print(kamodo.signatures)
        raise


def test_unit_function_composition():
    kamodo = Kamodo('X[m] = x', verbose=True)

    @kamodofy(units='km/s', arg_units=dict(x = 'm'))
    def v(x):
        return x

    kamodo['v'] = v
    kamodo['speed'] = 'v(X)'
    assert kamodo.speed.meta['units'] == 'km/s'


def test_method_args():
    class TestModel(Kamodo):
        def __init__(self, *args, **kwargs):
            super(TestModel, self).__init__(*args, **kwargs)
            self['rho'] = self.density

        def density(self, alt, lat, lon):
            return alt + lat + lon

        density.meta = dict(units='1/cm^3')

    test = TestModel(verbose=True)
    assert str(list(test.keys())[0].args[0]) != 'self'

    test['r(alt, lat, lon)[1/m^3]'] = 'rho'
    try:  # check that args are in the correct order
        assert list(test.signatures.values())[0]['symbol'].args == list(test.signatures.values())[1]['symbol'].args
    except:
        print('\n', test.detail())
        raise


def test_komodofy_decorator():
    @kamodofy(units='kg/cm^3')
    def my_density(x, y, z):
        return x + y + z

    assert my_density.meta['units'] == 'kg/cm^3'
    assert my_density(3, 4, 5) == 12


def test_komodofy_register():
    @kamodofy(units='kg/cm^3')
    def my_density(x, y, z):
        return x + y + z

    kamodo = Kamodo(rho=my_density)
    assert kamodo.rho.meta['units'] == 'kg/cm^3'


def test_komodofy_method():
    class TestModel(Kamodo):
        def __init__(self, *args, **kwargs):
            super(TestModel, self).__init__(*args, **kwargs)
            self['rho'] = self.density

        @kamodofy(units='1/cm^3')
        def density(self, alt, lat, lon):
            return alt + lat + lon

    test = TestModel(verbose=False)
    assert test.density.meta['units'] == '1/cm^3'


#   try:
#       assert len(kamodo.detail()) == 1
#   except:
#       print kamodo.symbol_registry
#       print kamodo.detail()
#       raise

# # def test_remove_symbol():
# #     kamodo = Kamodo(f = 'x')
# #     kamodo['g'] = '2*f'
# #     kamodo.remove_symbol('f')
# #     try:
# #         assert len(kamodo.symbol_registry) == 1
# #         assert len(kamodo) == 2
# #     except:
# #         print '\n',kamodo.detail()
# #         raise

# def test_function_arg_ordering():
#   def x(r,theta,phi):
#       '''converts from spherical to cartesian'''
#       return r*np.sin(theta)*np.cos(phi)

#   kamodo = Kamodo(x = x)
#   rhs_args = get_function_args(kamodo.x)
#   lhs_args = kamodo.signatures.values()[0]['symbol'].args

#   for i in range(3):
#       assert lhs_args[i] == rhs_args[i]

# def test_function_change_of_variables():
#   def rho(x):
#       return x**2

#   kamodo = Kamodo(x = 'y')
#   kamodo['rho'] = rho
#   kamodo['rho_perp'] = 'rho(x)'

#   assert str(get_function_args(kamodo.rho_perp)[0]) == 'y'


def test_kamodofy_update_attribute():
    @kamodofy(units='cm', update='y')
    def f(x):
        return x # pragma: no cover

    kamodo = Kamodo(f=f)
    assert f.update == 'y'


def test_kamodo_coupling():
    @kamodofy(units='cm', update='y')
    def y_iplus1(y, x):
        return y + x

    @kamodofy(units='m')
    def x_iplus1(x, y):
        return x - y

    kamodo = Kamodo()
    kamodo['x_iplus1'] = x_iplus1
    kamodo['y_iplus1'] = y_iplus1

    kamodo.x_iplus1.update = 'x'
    assert kamodo.x_iplus1.update == 'x'
    assert kamodo.y_iplus1.update == 'y'

    simulation = kamodo.simulate(y=1, x=0, steps=1)

    for state in simulation:
        pass

    assert state['x'] == -1
    assert state['y'] == 0


def test_units():
    nT = get_unit('nT')
    assert str(nT) == 'nanotesla'
    assert str(nT.abbrev) == 'nT'

    R_E = get_unit('R_E')
    assert str(R_E) == 'earth radii'
    assert str(R_E.abbrev) == 'R_E'


def test_default_composition():
    ## Need to wrap function defaults
    def create_wrapped(args, expr):
        func = lambdify(args, expr)

        @functools.wraps(func)
        def wrapped(*_args, **kwargs):
            # Write the logic here to parse _args and kwargs for the arguments as you want them
            return func(*actual_args)

        return wrapped


def test_vectorize():
    @np.vectorize
    @kamodofy(units='kg')
    def f(x, y):
        return x + y

    kamodo = Kamodo(f=f)
    kamodo.f([3, 4, 5], 5)


def test_redefine_variable():
    kamodo = Kamodo(rho='x + y')
    kamodo['rho'] = 'a + b'
    kamodo['rho(a,b)'] = 'a*b'

def test_unit_composition_registration():
    server = Kamodo(**{'M': kamodofy(lambda r=3: r, units='kg'),
                       'V[m^3]': (lambda r: r**3)}, verbose=True)
    user = Kamodo(mass=server.M, vol=server.V,
              **{'rho(r)[g/cm^3]':'mass/vol'}, verbose=True)

    result = (3/3**3)*(1000)*(1/10**6)
    assert np.isclose(user.rho(3), result)


def test_unit_expression_registration():
    kamodo = Kamodo(verbose=True)
    kamodo['f(x[cm])[cm**2]'] = 'x**2'


def test_multi_unit_composition():
    kamodo = Kamodo('a(x[s])[km] = x', verbose=True)
    kamodo['b(x[cm])[g]'] = 'x'
    kamodo['c'] = 'b(a)'
    print(kamodo.c.meta)
    assert kamodo.c.meta['units'] == 'g'
    assert kamodo.c.meta['arg_units']['x'] == str(get_abbrev(get_unit('s')))

def test_unit_composition_conversion():
    kamodo = Kamodo('a(x[kg])[m] = x', verbose=True)
    kamodo['b(x[cm])[g]'] = 'x'
    kamodo['c'] = 'b(a)'

    assert kamodo.c.meta['units'] == 'g'
    assert kamodo.c.meta['arg_units']['x'] == 'kg'
    assert kamodo.c(3) == kamodo.b(100*kamodo.a(3))


def test_get_arg_units():

    def assign_unit(symbol, **units):
        if isinstance(symbol, str):
            symbol = sympify(symbol)
        return symbol.subs(units)

    f_units = assign_unit('f(x,y)', x=get_unit('s'), y=get_unit('hour'))
    g_units = assign_unit('g(a,b,c)', a=get_unit('km'), b=get_unit('cm'), c=get_unit('m'))
    unit_registry = {
        sympify('f(x,y)'): f_units,
        f_units: get_unit('kg/m^3'),
        sympify('g(a,b,c)'): g_units,
        g_units: get_unit('m^2')
    }
    x, c = symbols('x c')
    result = get_arg_units(sympify('h(g(f(x,y)))'), unit_registry)
    assert result[x] == get_unit('s')
    result = get_arg_units(sympify('f(x,y)*g(a,b,c)'), unit_registry)
    assert result[c] == get_unit('m')

def test_compose_unit_multiply():
    kamodo = Kamodo('a(x[kg])[m] = x', verbose=True)
    kamodo['e'] = '2*a'



