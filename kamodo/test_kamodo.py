"""
Tests for kamodo.py

"""

import numpy as np
from sympy import symbols, Symbol
import pytest
from sympy.core.function import UndefinedFunction
import pandas as pd
from kamodo import Kamodo, get_unit, kamodofy, Eq
import functools
from sympy import lambdify, sympify
from kamodo import get_abbrev
from .util import get_arg_units
from .util import get_unit_quantity, convert_unit_to
from kamodo import from_kamodo, compose
from sympy import Function
from kamodo import KamodoAPI
from .util import serialize, NumpyArrayEncoder
from .util import get_kamodo_unit_system

import warnings


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

def test_Kamodo_reassignment_units():
    kamodo = Kamodo(verbose=True)
    kamodo['s(x[km],y[km])[kg]'] = 'x + y'
    assert kamodo.s(1, 1) == 2
    kamodo['s(x[m],y[m])[g]'] = '3*x + y'
    assert kamodo.s(1, 1) == 4


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
    assert 'h' in kamodo.signatures


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
    warnings.simplefilter('error')
    kamodo = Kamodo(f='x**2', verbose=True)
    assert str(kamodo.to_latex(mode='inline')) == r'$f{\left(x \right)} = x^{2}$'
    kamodo = Kamodo(g='x', verbose=True)
    assert str(kamodo.to_latex(mode='inline')) == r'$g{\left(x \right)} = x$'
    kamodo['f(x[cm])[kg]'] = 'x**2'
    kamodo['g'] = kamodofy(lambda x: x**2, units='kg', arg_units=dict(x='cm'), equation='$x^2$')
    kamodo['h'] = kamodofy(lambda x: x**2, units='kg', arg_units=dict(x='cm'))
    
    @kamodofy(units = 'kg/m^3', citation = 'Bob et. al, 2018')
    def rho(x = np.array([3,4,5]), y = np.array([1,2,3])):
        """A function that computes density"""
        return x+y
    kamodo['rho'] = rho
    kamodo.to_latex()


def test_expr_conversion():
    kamodo = Kamodo('$a[km] = x$', verbose=True)
    print(kamodo.items())
    kamodo.a

def test_get_unit_fail():
    with pytest.raises(NameError):
        get_unit('unregistered units')
    with pytest.raises(NameError):
        get_unit('runregistered')

def test_get_unit_quantity():
    mykm = get_unit_quantity('mykm', 'km', scale_factor=2)
    mygm = get_unit_quantity('mygm', 'gram', scale_factor=4)
    assert str(convert_unit_to(mykm, get_unit('m'))) == '2000*meter'
    assert str(convert_unit_to(mygm, get_unit('kg'))) == 'kilogram/250'

# def test_validate_units():
#     f, x = symbols('f x')
#     with pytest.raises(ValueError):
#         validate_units(Eq(f, x), dict(f=get_unit('kg'), x=get_unit('m')))
#     validate_units(Eq(f, x), dict(f=get_unit('kg'), x=get_unit('g')))

#     lhs_units = validate_units(sympify('f(x)'), dict(f=get_unit('kg'), x=get_unit('m')))
#     print(lhs_units)


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
        assert get_unit(kamodo.signatures['p']['units']) == get_unit('kg*km/s')
    except:
        print(kamodo.signatures)
        raise

def test_unit_composition_mixed():
    kamodo = Kamodo('$rho[kg/m^3] = x^3$', '$v[cm/s] = y^2$', verbose=True)
    kamodo['p[Pa]'] = '$\\rho v^2$'


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
    kamodo.evaluate('f', x=3,y=2)


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


def test_compose_unit_add():
    kamodo = Kamodo(verbose=True)

    @kamodofy(units='m', arg_units={'x': 'kg'})
    def a(x):
        return x

    kamodo['a'] = a
    kamodo['b(y[cm])[km]'] = 'y'
    kamodo['c(x,y)[km]'] = '2*a + 3*b'
    assert kamodo.c.meta['arg_units']['x'] == str(get_abbrev(get_unit('kg')))
    assert kamodo.c.meta['arg_units']['y'] == str(get_abbrev(get_unit('cm')))
    result = 2*(3)/1000 + 3*(3)
    assert kamodo.c(3, 3) == result

def test_compose_unit_raises():

    with pytest.raises(NameError):
        kamodo = Kamodo('a(x[kg])[m] = x',
                        'b(y[cm])[km] = y', verbose=True)

        # this should fail since x is registered with kg
        kamodo['c(x[cm],y[cm])[km]'] = '2*a + 3*b'


def test_repr_latex():
    kamodo = Kamodo(f='x')
    assert kamodo._repr_latex_() == r'\begin{equation}f{\left(x \right)} = x\end{equation}'
    kamodo = Kamodo(f=lambda x:x)
    assert kamodo.f._repr_latex_() == '$f{\\left(x \\right)} = \\lambda{\\left(x \\right)}$'


def test_dataframe_detail():
    kamodo = Kamodo(f='x')
    type(kamodo.detail())
    assert len(kamodo.detail()) == 1
    assert isinstance(kamodo.detail(), pd.DataFrame)

def test_from_kamodo():
    kamodo = Kamodo(f='x')
    knew = from_kamodo(kamodo, g='f**2')
    assert knew.g(3) == 9


def test_jit_evaluate():
    """Just-in-time evaluation"""
    kamodo = Kamodo(f='x')
    assert kamodo.evaluate('g=f**2', x = 3)['x'] == 3
    with pytest.raises(SyntaxError):
        kamodo.evaluate('f**2', x = 3)['x'] == 3

def test_eval_no_defaults():
    kamodo = Kamodo(f='x', verbose=True)
    kamodo['g'] = lambda x=3: x
    kamodo['h'] = kamodofy(lambda x=[2,3,4]: (kamodofy(lambda x_=np.array([1]): x_**2) for x_ in x))
    assert kamodo.evaluate('f', x=3)['x'] == 3
    assert kamodo.evaluate('g')['x'] == 3
    assert kamodo.evaluate('h')['x_'][-2] == 1.

def test_compose():
    k1 = Kamodo(f='x')
    k2 = Kamodo(g='y**2', h = kamodofy(lambda x: x**3))
    k3 = compose(m1=k1, m2=k2)
    assert k3.f_m1(3) == 3
    assert k3.g_m2(3) == 3**2
    assert k3.h_m2(3) == 3**3
    k3['h(f_m1)'] = 'f_m1'
    assert k3.h(3) == 3

def test_symbol_replace():
    k = Kamodo(f='x', verbose=True)

    f1, f2 = list(k.keys())
    print('\n|||||',f1, f2, '||||||')
    k[f1] = 'x**2'
    assert k.f(3) == 9
    print('\n|||||', *list(k.keys()), '|||||')
    k[f2] = 'x**3'
    print('\n|||||', *list(k.keys()), '|||||')
    assert k.f(3) == 27

def test_contains():
    kamodo = Kamodo(f='x', verbose=True)
    assert 'f' in kamodo
    assert 'f( x )' in kamodo
    f, x = symbols('f x')
    assert f in kamodo # f is a symbo
    f = Function(f) # f is an Unefined Function
    assert f in kamodo
    assert f(x) in kamodo
    assert f('x') in kamodo

def test_unusual_signature():
    with pytest.raises(NotImplementedError):
        kamodo = Kamodo()
        kamodo['f(x)=f(cm)=kg'] = 'x'
    with pytest.raises(NotImplementedError):
        kamodo = Kamodo('f(x)=f(cm)=kg=x')


# def test_remove_symbol():
#     kamodo = Kamodo(f='x', verbose=True)
#     kamodo.remove_symbol('f')
#     assert 'f' not in kamodo

def test_method_registry():
    
    class MyClass(Kamodo):
        @kamodofy
        def f(self, x):
            return x**2

    myclass = MyClass()
    myclass['f'] = myclass.f

def test_del_function():
    kamodo = Kamodo(f='x', g='y', h='y', verbose=True)
    del(kamodo.f)
    assert 'f' not in kamodo
    assert 'f' not in kamodo.signatures
    del(kamodo['g'])
    assert 'g' not in kamodo
    del(kamodo['h(y)'])
    print(kamodo.keys())
    assert 'h(y)' not in kamodo

    with pytest.raises(AttributeError):
        del(kamodo.y)

def test_simple_figure():
    @kamodofy(units='kg', hidden_args=['ions'])
    def f_N(x_N):
        return x_N**2

    kamodo = Kamodo(f_N=f_N,verbose=True)
    kamodo.plot(f_N=dict(x_N=np.linspace(-4, 3, 30)))

def test_unavailable_4d_plot_type():
    def g(x=np.array([1]),
          y=np.array([1]),
          z=np.array([1]),
          t=np.array([1])):
        return x**2 + y**2 + z**2 + t**2

    kamodo = Kamodo(g=g, verbose=True)
    with pytest.raises(KeyError):
        kamodo.plot('g')

def test_multiple_traces():
    kamodo = Kamodo(f='x', g='x**2')
    kamodo.plot(
        f=dict(x=np.linspace(-1, 1, 10)),
        g=dict(x=np.linspace(-5, 5, 10)))

def test_unitless_composition():
    @kamodofy
    def alpha(x):
        return x

    @kamodofy
    def beta(y):
        return y


    kamodo = Kamodo(alpha=alpha, beta_=beta, verbose=True)
    kamodo['Gamma'] = 'alpha(beta_)'
    kamodo

def test_reserved_name():
    kamodo = Kamodo(verbose=True)
  
    @kamodofy
    def test(x, y):
        return x+y
    with pytest.raises(NotImplementedError):
        kamodo['test'] = test


class Ktest(Kamodo):
    def __init__(self, **kwargs):
        super(Ktest, self).__init__()

        t_N = pd.date_range('Nov 9, 2018', 'Nov 20, 2018', freq = 'H')
        @kamodofy(units='kg/m^3')
        def rho_N(t_N=t_N):
            t_N = pd.DatetimeIndex(t_N)
            t_0 = pd.to_datetime('Nov 9, 2018') 
            try:
                dt_days = (t_N - t_0).total_seconds()/(24*3600)
            except TypeError as err_msg:
                return 'cannot work with {} {}  {}'.format(type(t_N), type(t_N[0]), err_msg)

            result = np.abs(weierstrass(dt_days))
            return result

        @kamodofy(units='nPa')
        def p(x=np.linspace(-5, 5, 30)):
            try:
                return x**2
            except TypeError as m:
                print(m)
                print(type(x), x[0])
                raise

        @kamodofy(
            equation=r"\sum_{n=0}^{500} (1/2)^n cos(3^n \pi x)",
            citation='https://en.wikipedia.org/wiki/Weierstrass_function'
            )
        def weierstrass(x = np.linspace(-2, 2, 1000)):
            '''
            Weierstrass  function
            A continuous non-differentiable 
            https://en.wikipedia.org/wiki/Weierstrass_function
            '''
            nmax = 500
            n = np.arange(nmax)

            xx, nn = np.meshgrid(x, n)
            ww = (.5)**nn * np.cos(3**nn*np.pi*xx)
            return ww.sum(axis=0)
                

        self['rho_N'] = rho_N
        self['p'] = p
        self['Weierstrass'] = weierstrass


def test_kamodo_inline_merge():
    k1 = Kamodo(f='x**2')
    k2 = Kamodo(g=lambda y: y-1)

    # create a local namespace holding both kamodo objects
    ns = {'k1':k1, 'k2': k2}
    k3 = Kamodo(myf = sympify('k1.f(x) + k2.g(y)', locals=ns))
    assert k3.myf(x=3, y=4) == 3**2 + 4 - 1

def test_default_forwarding():
    x = np.linspace(-5, 5, 12)
    
    def f(x=x):
        return x**2
    
    k = Kamodo(f=f)
    k['g'] = 'f+2'
    assert len(k.g()) == 12
    
def test_multi_arg_units():
    kamodo = Kamodo(verbose=True)

    # @kamodofy(units='m', arg_units={'X': 'kg', 'Y':'cm', 'Z': 's'})
    # def f(X, Y, Z):
    #     return x*y*z

    kamodo['f(X[kg],Y[cm],Z[s])[m]'] = 'X*Y*Z'
    kamodo['a[g]'] = 'x'
    kamodo['b[m]'] = 'y'
    kamodo['c[ms]'] = 'z**2'
    kamodo['d(x,y,z)[cm]'] = 'f(a,b,c)'

def test_broken_unit():
    k = Kamodo()
    k['f[N]'] = 'x'

    get_unit('newton')
    get_unit('N')

def test_frequency_composition():
    @kamodofy(units='rad/s', arg_units={'B':'T', 'n_e':'1/m**3'})
    def omega_uh1(B, n_e):
        return np.sqrt(B**2+n_e**2)


    kamodo_test = Kamodo(verbose=True)
    kamodo_test['B_mag'] = kamodofy(lambda B=np.linspace(0.1,1.,10): B, units='nT', arg_units={'B':'nT'})
    kamodo_test['n_e'] = kamodofy(lambda n=np.linspace(4.,13.,10)*10**19:n, units='1/m**3', arg_units={'n':'1/m**3'})
    kamodo_test['omega_uh1'] = omega_uh1
    print(kamodo_test.unit_registry)

    #---------(input)--------
    kamodo_test['omega_uh1A'] = 'omega_uh1(B_mag, n_e)'
    kamodo_test.omega_uh1A


def test_frequency_units():
    omega = get_unit('rad')/get_unit('s')
    freq = get_unit('deg')/get_unit('s')
    kamodo_units = get_kamodo_unit_system()
    convert_unit_to(omega, freq, kamodo_units)

