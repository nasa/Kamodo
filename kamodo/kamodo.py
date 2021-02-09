"""
Copyright © 2017 United States Government as represented by the Administrator, National Aeronautics and Space Administration.  
No Copyright is claimed in the United States under Title 17, U.S. Code.  All Other Rights Reserved.
"""

import numpy as np
from sympy import Integral, Symbol, symbols, Function

try:
    from sympy.parsing.sympy_parser import parse_expr
except ImportError:  # may occur for certain versions of sympy
    from sympy import sympify as parse_expr

from collections import OrderedDict
from collections import UserDict
import collections

from sympy import lambdify
from sympy.parsing.latex import parse_latex
from sympy import latex
from sympy.core.function import UndefinedFunction
from inspect import getfullargspec
from sympy import Eq
import pandas as pd
# from IPython.display import Latex


from sympy.physics import units as sympy_units
from sympy.physics.units import Quantity
from sympy.physics.units import Dimension
from sympy import Expr

import functools
from .util import kamodofy
from .util import sort_symbols
from .util import simulate
from .util import unit_subs
from .util import get_defaults, valid_args, eval_func
# from .util import to_arrays, cast_0_dim
from .util import beautify_latex, arg_to_latex
from .util import concat_solution
from .util import convert_to
from .util import unify, get_abbrev, get_expr_unit
from .util import is_function, get_arg_units

import sympy.physics.units as u

import plotly.graph_objs as go
from plotly import figure_factory as ff

from plotting import plot_dict, get_arg_shapes, get_plot_key
from .util import existing_plot_types
from .util import get_dimensions
from .util import reserved_names

from sympy import Wild
from types import GeneratorType
import inspect

import re

import urllib.request, json
import requests
from .util import serialize, deserialize
from .util import sign_defaults
import forge



def get_unit_quantities():
    """returns a map of sympy units {symbol: quantity} """
    subs = {}
    for k, v in list(sympy_units.__dict__.items()):
        if isinstance(v, Expr) and v.has(sympy_units.Unit):
            subs[Symbol(k)] = v  # Map the `Symbol` for a unit to the quantity
    return subs


def clean_unit(unit_str):
    """remove brackets and all white spaces in and around unit string"""
    return unit_str.strip().strip('[]').strip()


def get_unit(unit_str, unit_subs=unit_subs):
    """get the unit quantity corresponding to this string

    unit_subs should contain a dictonary {symbol: quantity} of custom units not available in sympy"""

    unit_str = clean_unit(unit_str)
    units = get_unit_quantities()

    if unit_subs is not None:
        units.update(unit_subs)

    if len(unit_str) == 0:
        return Dimension(1)
    try:
        unit = parse_expr(unit_str.replace('^', '**')).subs(units)
    except:
        raise NameError('something wrong with unit str [{}], type {}'.format(unit_str, type(unit_str)))
    try:
        assert len(unit.free_symbols) == 0
    except:
        raise NameError("Unsupported unit: {} {} {}".format(
            unit_str, type(unit_str), unit.free_symbols))
    return unit




def args_from_dict(expr, local_dict, verbose):
    """retrieve the symbols of an expression

    if a  {str: symbol} mapping is provided, use its symbols instead
    """
    args = []
    # if local_dict is not None:
    if verbose:
        print('args_from_dict: available names {}'.format(local_dict.keys()))
    for a in expr.args:
        if verbose:
            print('args_from_dict: {}'.format(a))
        if str(a) in local_dict:
            args.append(local_dict[str(a)])
        else:
            args.append(a)
    return tuple(args)
    # else:
    #     return expr.args

def get_str_units(bracketed_str):
    """finds all bracketed units in a string

    supports functions like
        bracketed_str = 'f(x[m],y[cm**3],z[km^.5],q[(cm*g)^-3])[km/s]'
        get_str_units(bracketed_str)
        >>> ['m', 'cm**3', 'km^.5', '(cm*g)^-3', 'km/s']
    """
    return re.findall(r"\[([A-Za-z0-9_/*^\\(\\)-\.]+)\]", bracketed_str)

def str_has_arguments(func_str): 
    bracket = func_str.find('[')
    paren = func_str.find('(')
    if bracket > 0:
        if paren > 0:
            # f(x[(cm)^2]) has arguments since paren comes first
            # f[(cm)^2] has no arguments since bracket comes first
            return bracket > paren
        return False # no paren hence no arguments
    # no bracket
    if paren > 0:
        # has argument
        return True
    return False


def extract_units(func_str):
    """separate the units from the left-hand-side of an expression assuming bracket notation

    We return symbol and dictionary mapping symbols to units

    needs to handle the following cases:

        args and output have units
            extract_units('f(x[cm],y[km])[kg]')
            ('f(x,y)', {'x': 'cm', 'y': 'km', 'f(x,y)': 'kg'})

        args have parenthesis in units and output has no units
            extract_units('f(x[(cm )^2])')
            ('f(x)', {'x': '(cm)^2', 'f(x)': ''})

        output has parenthesis in units
            extract_units('f[(cm)^2]')
            ('f', {'f': '(cm)^2'})

        no args named and no units named
            ('f', {'f': ''})

    """
    # remove any spaces
    func_str = func_str.replace(' ', '')
    all_units = get_str_units(func_str)

    unit_dict = {}
    if str_has_arguments(func_str):
        # f(x[cm],y[km])[kg] -> x[cm],y[km]
        lhs_args = re.findall(r"\((.+)\)", func_str)[0]

        # 'x[cm],y[km]' -> ['cm', 'km']
        arg_units = get_str_units(lhs_args)
        args = lhs_args.split(',')

        for arg, unit in zip(args, arg_units):
            unit_dict[arg.replace('[{}]'.format(unit), '')] = unit

        if len(all_units) == len(arg_units):
            output_units = ''
        else:
            output_units = all_units[-1]
    else:
        if len(all_units) > 0:
            output_units = all_units[0]
        else:
            output_units = ''


    # f(x[cm],y[km])[kg]
    lhs = func_str
    for unit in all_units:
        lhs = lhs.replace('[{}]'.format(unit), '')

    unit_dict[lhs] = output_units
    return lhs, unit_dict


def expr_to_symbol(expr, args):
    """f,args -> f(args) of class function, f(x) -> f(x)"""
    if type(expr) == Symbol:
        return parse_expr(str(expr) + str(args))
        # return symbols(str(expr), cls = Function)
    else:
        return expr


def parse_lhs(lhs, local_dict, verbose):
    """parse the left-hand-side

    returns: symbol, arguments, units, parsed_expr
    """
    lhs, unit_dict = extract_units(lhs)
    if lhs in reserved_names:
        raise NameError('{} is a reserved name'.format(lhs))
    parsed = parse_expr(lhs)
    args = args_from_dict(parsed, local_dict, verbose)
    symbol = expr_to_symbol(parsed, args)
    return symbol, args, unit_dict, parsed


def parse_rhs(rhs, is_latex, local_dict):
    """parse the right-hand-side of an equation

    subsititute variables from local_dict where available
    """
    if is_latex:
        expr = parse_latex(rhs).subs(local_dict)
    else:
        try:
            expr = parse_expr(rhs).subs(local_dict)
        except SyntaxError:
            print('cannot parse {} with {}'.format(rhs, local_dict))
            raise
    return expr


def get_function_args(func, hidden_args=[]):
    """converts function arguments to list of symbols"""
    return symbols([a for a in getfullargspec(func).args if a not in hidden_args])




# class Kamodo(collections.OrderedDict):
class Kamodo(UserDict):
    '''Kamodo base class demonstrating common API for space weather models


    This API provides access to space weather fields and their properties through:

        interpolation of variables at user-defined points
        unit conversions
        coordinate transformations specific to space weather domains

    Required methods that have not been implemented in child classes
    will raise a NotImplementedError
    '''

    def __init__(self, *funcs, **kwargs):
        """Base initialization method

        Args:
            param1 (str, optional): Filename of datafile to interpolate from

        """

        super(Kamodo, self).__init__()
        self.symbol_registry = OrderedDict()
        self.unit_registry = OrderedDict()

        symbol_dict = kwargs.pop('symbol_dict', None)

        self.verbose = kwargs.pop('verbose', False)
        self.signatures = OrderedDict()

        for func in funcs:
            if self.verbose:
                print('registering {}'.format(func))
            if type(func) is str:
                components = func.split('=')
                if len(components) == 2:
                    # function(arg)[unit] = expr
                    lhs, rhs = components
                    self[lhs.strip('$')] = rhs
                else:
                    raise NotImplementedError(
                        'cannot register functions of the form {}'.format(func))

        for sym_name, expr in list(kwargs.items()):
            if self.verbose:
                print('registering {} with {}'.format(sym_name, expr))
            self[sym_name] = expr


    def register_symbol(self, symbol):
        self.symbol_registry[str(type(symbol))] = symbol


    # def remove_symbol(self, sym_name):
    #     """remove the sym_name (str) from the object"""

    #     if self.verbose:
    #         print('remove_symbol: removing {} from symbol_registry'.format(sym_name))
    #     try:
    #         symbol = self.symbol_registry.pop(sym_name)
    #     except KeyError:
    #         raise KeyError('{} was not in symbol_registry {}'.format(
    #             sym_name, self.symbol_registry.keys()))
    #     if self.verbose:
    #         print('remove_symbol: removing {} from signatures'.format(symbol))
    #     self.signatures.pop(str(type(symbol)))
    #     if self.verbose:
    #         print('remove_symbol: removing {} from unit_registry'.format(symbol))
    #     remove = []
    #     for sym in self.unit_registry:
    #         if is_function(sym): # {rho(x): rho(cm), rho(cm): kg}
    #             if type(sym) == type(symbol):
    #                 remove.append(sym)
    #     for sym in remove:
    #         self.unit_registry.pop(sym)

    #     if self.verbose:
    #         print('removing {} from {}'.format(symbol, self.keys()))

    #     if symbol in self.keys():
    #         try:
    #             self.pop(symbol)
    #         except KeyError:
    #             if self.verbose:
    #                 print('something wrong with {} in {}'.format(symbol, self.keys()))
    #             raise KeyError('{} not in {}'.format(symbol, self.keys()))

    #     if type(symbol) in self.keys():
    #         try:
    #             self.pop(type(symbol))
    #         except KeyError:
    #             raise KeyError('{} not in {}'.format(type(symbol), self.keys()))


    def parse_key(self, sym_name):
        """parses the symbol name

        sym_name must be a string

        returns: symbol, args, unit_dict, lhs_expr
        """
        args = tuple()
        try:
            symbol, args, unit_dict, lhs_expr = parse_lhs(
                sym_name,
                self.symbol_registry,
                self.verbose)
        except:
            raise NotImplementedError('could not parse key {}'.format(sym_name))

        # try:
        if len(args) > 0:
            symbol_str = str(symbol).replace(' ', '')
        else:
            symbol_str = str(type(symbol))
        units = get_abbrev(unit_dict[symbol_str])
        # except KeyError:
        #     raise KeyError('{} not in {}'.format(symbol_str, unit_dict))

        return symbol, args, units, lhs_expr

    def parse_value(self, rhs_expr, local_dict):
        """returns an expression from string"""
        if type(rhs_expr) is str:
            is_latex = '$' in rhs_expr
            rhs_expr = rhs_expr.strip('$').strip()
            rhs_expr = parse_rhs(rhs_expr, is_latex, local_dict)
        return rhs_expr

    def check_or_replace_symbol(self, symbol, free_symbols, rhs_expr=None):
        """Rules of replacement:

        """
        if self.verbose:
            print('symbol arguments: >>> {} {} <<<'.format(symbol.args, type(symbol.args)))
            print('free_symbols: >>> {} {} <<<'.format(free_symbols, type(free_symbols)))
        # try:
        lhs_args = set(symbol.args)
        # except TypeError:
        #     lhs_args = free_symbols
        #     symbol = Function(str(symbol))(*lhs_args)

        if lhs_args != set(free_symbols):
            free_symbols_ = tuple(free_symbols)
            if len(free_symbols) == 1:
                # try:
                symbol = parse_expr(str(type(symbol)) + str(free_symbols_))
                # except:
                #     raise NotImplementedError('cannot parse', str(type(symbol)) + str(free_symbols_))
            else:
                if len(lhs_args) > 0:
                    raise NameError("Mismatched symbols {} and {}".format(lhs_args, free_symbols))
                else:
                    if self.verbose:
                        print('type of rhs symbols:', type(free_symbols))
                    # if isinstance(free_symbols, set):
                    #     free_symbols_ = sort_symbols(free_symbols)
                    if self.verbose:
                        print('replacing {} with {}'.format(
                            symbol, str(type(symbol)) + str(free_symbols_)))
                    # try:
                    symbol = parse_expr(str(type(symbol)) + str(free_symbols_))
                    # except:
                    #     raise NotImplementedError('cannot parse', str(type(symbol)) + str(free_symbols_))

        return symbol

    def validate_function(self, lhs_expr, rhs_expr):
        assert lhs_expr.free_symbols == rhs_expr.free_symbols


    def vectorize_function(self, symbol, rhs_expr, composition):
        try:
            func = lambdify(symbol.args, rhs_expr, modules=['numexpr'])
            if self.verbose:
                print('lambda {} = {} labmdified with numexpr'.format(symbol.args, rhs_expr))
        except: # numexpr not installed
            func = lambdify(symbol.args, rhs_expr, modules=['numpy', composition])
            if self.verbose:
                print('lambda {} = {} lambdified with numpy and composition:'.format(
                    symbol.args, rhs_expr))
                for k, v in composition.items():
                    print('\t', k, v)
        signature = sign_defaults(symbol, rhs_expr, composition)
        return signature(func)


    def update_unit_registry(self, func, arg_units):
        """Inserts unit functions into registry"""
        lhs, unit_dict = extract_units(func)
        # if arg_units is None:
        #     arg_units = {}
        for key, value in unit_dict.items():
            if key != lhs:
                arg_units[parse_expr(key)] = get_unit(value)

        lhs_expr = parse_expr(lhs)
        func_units = lhs_expr.subs(arg_units)
        self.unit_registry[lhs_expr] = func_units
        output_units = unit_dict[lhs]
        self.unit_registry[func_units] = get_unit(output_units)
        return lhs


    def register_signature(self, symbol, units, lhs_expr, rhs_expr):
        # if isinstance(units, str):
        unit_str = units
        if self.verbose:
            print('unit str {}'.format(unit_str))

        self.signatures[str(type(symbol))] = dict(
            symbol=symbol,
            units=unit_str,
            lhs=lhs_expr,
            rhs=rhs_expr,
            # update = getattr(self[lhs_expr],'update', None),
        )

    def register_function(self, func, lhs_symbol, lhs_expr, lhs_units):
        hidden_args = []
        if hasattr(func, 'meta'):
            hidden_args = func.meta.get('hidden_args', [])

        if type(func) is np.vectorize:
            rhs_args = get_function_args(func.pyfunc, hidden_args)
        else:
            rhs_args = get_function_args(func, hidden_args)
        if str(rhs_args[0]) == 'self':  # in case function is a class method
            rhs_args.pop(0)

        lhs_symbol_old = lhs_symbol
        lhs_symbol = self.check_or_replace_symbol(lhs_symbol, rhs_args)
        units = lhs_units
        if hasattr(func, 'meta'):
            if self.verbose:
                print('function has metadata', func.meta.keys())
            if len(lhs_units) > 0:
                if lhs_units != func.meta['units']:
                    raise NameError('Units mismatch:{} != {}'.format(lhs_units, func.meta['units']))
            else:
                units = func.meta['units']
            rhs = func.meta.get('equation', func)
            if self.verbose:
                print('rhs from meta: {}'.format(rhs))
        else:
            rhs = func
            if self.verbose:
                print('rhs from input func {}'.format(rhs))
            try:
                setattr(func, 'meta', dict(units=lhs_units, arg_units=None))
            except:  # will not work on bound methods
                pass

        arg_units = func.meta.get('arg_units', None)
        if arg_units is not None:
            unit_dict = {arg: get_unit(arg_unit) for arg, arg_unit in arg_units.items()}
            unit_expr = lhs_symbol.subs(unit_dict)
            self.unit_registry[lhs_symbol] = unit_expr
            self.unit_registry[unit_expr] = get_unit(units)
        else:
            self.unit_registry[lhs_symbol] = get_unit(units)

        self.register_signature(lhs_symbol, units, lhs_expr, rhs)
        try:
            func._repr_latex_ = lambda: self.func_latex(str(type(lhs_symbol)), mode='inline')
        except AttributeError:
            # happens with bound methods
            pass
        super(Kamodo, self).__setitem__(lhs_symbol, func)  # assign key 'f(x)'
        super(Kamodo, self).__setitem__(type(lhs_symbol), self[lhs_symbol])  # assign key 'f'
        self.register_symbol(lhs_symbol)


    def __setitem__(self, sym_name, input_expr):
        """Assigns a function or expression to a new symbol,
        performs unit conversion where appropriate

        """
        if not isinstance(sym_name, str):
            sym_name = str(sym_name)

        symbol, args, lhs_units, lhs_expr = self.parse_key(sym_name)

        if hasattr(input_expr, '__call__'):
            self.register_function(input_expr, symbol, lhs_expr, lhs_units)

        else:
            if self.verbose:
                print(
                    "\n\nPARSING WITH UNIFY",
                    lhs_expr,
                    symbol,
                    lhs_units,
                    len(lhs_units),
                    type(lhs_units))
                print('symbol registry:', self.symbol_registry)

            rhs_expr = self.parse_value(input_expr, self.symbol_registry)

            if self.verbose:
                print('parsed rhs_expr', rhs_expr)

            if not isinstance(symbol, Symbol):
                if isinstance(lhs_expr, Symbol):
                    symbol = Function(lhs_expr)(*tuple(rhs_expr.free_symbols))
                else: #lhs is already a function
                    symbol = lhs_expr
                lhs_str = str(symbol)
                sym_name = sym_name.replace(str(lhs_expr), lhs_str)
            if self.verbose:
                print('unit registry contents:')
                for k, v in self.unit_registry.items():
                    print('\t', k, type(k), v)
            if '[' in sym_name:
                if self.verbose:
                    print('updating unit registry with {} -> {}'.format(sym_name, rhs_expr))
                rhs = rhs_expr
                arg_units = get_arg_units(rhs_expr, self.unit_registry)
                if self.verbose:
                    print(arg_units)
                sym_name = self.update_unit_registry(sym_name, arg_units)
                if self.verbose:
                    print('unit registry update returned', sym_name, self.unit_registry.get(symbol))
            else:

                if self.verbose:
                    print(sym_name,
                          symbol,
                          'had no units. Getting units from {}'.format(rhs_expr))

                expr_unit = get_expr_unit(rhs_expr, self.unit_registry, self.verbose)
                arg_units = get_arg_units(rhs_expr, self.unit_registry)


                if self.verbose:
                    print('registering {} with {} {}'.format(symbol, expr_unit, arg_units))

                if (symbol not in self.unit_registry) and (expr_unit is not None):
                    self.unit_registry[symbol] = symbol.subs(arg_units)
                    self.unit_registry[symbol.subs(arg_units)] = expr_unit


                if expr_unit is not None:
                    expr_dimensions = get_dimensions(expr_unit)
                    if expr_dimensions != Dimension(1):
                        lhs_units = str(get_abbrev(get_expr_unit(
                            expr_unit,
                            self.unit_registry,
                            self.verbose)))
                    else:
                        lhs_units = ''

                if self.verbose:
                    print('registered lhs_units', lhs_units)

                rhs = rhs_expr
                sym_name = str(sym_name)



            if len(lhs_units) > 0:
                if self.verbose:
                    print('about to unify lhs_units {} {} with {}'.format(
                        lhs_units, type(lhs_units), rhs))

                expr = unify(
                    Eq(parse_expr(sym_name), rhs),
                    self.unit_registry,
                    # to_symbol = symbol,
                    verbose=self.verbose)
                rhs_expr = expr.rhs

            if self.verbose:
                print('symbol after unify', symbol, type(symbol), rhs_expr)
                print('unit registry to resolve units:', self.unit_registry)

            units = get_expr_unit(symbol, self.unit_registry)
            if get_dimensions(units) != Dimension(1):
                units = get_abbrev(units)
                if units is not None:
                    units = str(units)
                else:
                    units = ''
            else:
                units = ''

            if self.verbose:
                print('units after resolve', symbol, units)
                for k, v in self.unit_registry.items():
                    print('\t{}: {}'.format(k, v))

            rhs_args = rhs_expr.free_symbols

            symbol = self.check_or_replace_symbol(symbol, rhs_args, rhs_expr)
            self.validate_function(symbol, rhs_expr)

            composition = {str(k_): self[k_] for k_ in self}
            arg_units = {}
            if symbol in self.unit_registry:
                unit_args = self.unit_registry[symbol]
                if unit_args is not None:
                    if len(unit_args.args) == len(symbol.args):
                        for arg, unit in zip(symbol.args, unit_args.args):
                            arg_units[str(arg)] = str(get_abbrev(unit))
            func = self.vectorize_function(symbol, rhs_expr, composition)
            meta = dict(units=units, arg_units=arg_units)
            func.meta = meta
            func.data = None
            self.register_signature(symbol, units, lhs_expr, rhs_expr)
            func._repr_latex_ = lambda: self.func_latex(str(type(symbol)), mode='inline')
            super(Kamodo, self).__setitem__(symbol, func)
            super(Kamodo, self).__setitem__(type(symbol), self[symbol])
            self.register_symbol(symbol)


    def __getitem__(self, key):
        try:
            return super(Kamodo, self).__getitem__(key)
        except KeyError:
            try:
                symbol = self.symbol_registry[str(key)]
            except KeyError:
                symbol = self.symbol_registry[str(type(key))]
            if symbol in self:
                return super().__getitem__(symbol)

    def __contains__(self, item):
        func_str = str(item).replace(' ', '')
        for key in self.keys():
            if func_str == str(key).replace(' ', ''):
                return True
        return False

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    def __delattr__(self, name):

        if name in self:
            symbol = self.symbol_registry[name]
            if self.verbose:
                print('__delattr__: removing name {}, symbol {}'.format(name, symbol))
                print(self.keys())
            # del self[symbol]
            super().pop(symbol)
            super().pop(type(symbol))

        else:
            raise AttributeError("No such field: " + name)

    def __delitem__(self, key):
        if self.verbose:
            print('__delitem__: got {} type: {}'.format(key, type(key)))
        if isinstance(key, str):
            try:
                symbol = self.symbol_registry[key]
            except KeyError:
                symbol, args, unit_dict, parsed = parse_lhs(
                    key, self.symbol_registry, self.verbose)
        else:
            symbol = key
        if self.verbose:
            print('__delitem__: removing {} {}'.format(symbol, type(symbol)))

        remove_keys = []
        for k in self.data:
            if str(k) == str(symbol):
                remove_keys.append(k)
            if str(k) == str(type(symbol)):
                remove_keys.append(k)
        if self.verbose:
            print('__delitem__: removing keys:', remove_keys)

        for key_ in remove_keys:
            self.data.pop(key_)

    # def get_units_map(self):
    #     """Maps from string units to symbolic units"""
    #     d = dict()
    #     star = Wild('star')
    #     for k, v in list(self.signatures.items()):
    #         unit = get_unit(v['units'])
    #         d[k] = unit
    #         d[v['symbol']] = unit

    #     return d

    def func_latex(self, key, mode='equation'):
        repr_latex = ""
        lhs = self.signatures[key]['symbol']
        rhs = self.signatures[key]['rhs']
        units = self.signatures[key]['units']
        arg_units = get_arg_units(lhs, self.unit_registry)

        if len(units) > 0:
            units = '{}'.format(get_abbrev(units))
        else:
            units = ''

        if len(arg_units) > 0:
            arg_strs = []
            for arg, arg_unit in arg_units.items():
                arg_strs.append("{}[{}]".format(latex(arg), get_abbrev(arg_unit)))
            lhs_str = "{}({})".format(latex(type(lhs)), ",".join(arg_strs))
        else:
            lhs_str = latex(lhs)
            # lhs_str = "{}({})".format(
            #     latex(type(lhs)),
            #     ','.join([latex(s) for s in lhs.args]))

        if len(units) > 0:
            lhs_str += "[{}]".format(units)

        latex_eq = ''
        latex_eq_rhs = ''

        if isinstance(rhs, str):
            latex_eq_rhs = rhs       
        elif hasattr(rhs, '__call__') | (rhs is None):
            lambda_ = symbols('lambda', cls=UndefinedFunction)
            # latex_eq = latex(Eq(lhs, lambda_(*lhs.args)), mode=mode)
            latex_eq_rhs = latex(lambda_(*lhs.args)) # no $$
        else:
            # latex_eq = latex(Eq(lhs, rhs), mode=mode)
            latex_eq_rhs = latex(rhs) # no $$ delimiter

        if len(str(units)) > 0:
            latex_eq = latex_eq.replace('=', '[{}] ='.format(units))

        if mode == 'equation':
            repr_latex += r"\begin{equation}"
            repr_latex += "{} = {}".format(lhs_str, latex_eq_rhs)
            repr_latex += r"\end{equation}"
        else:
            repr_latex += r"$"
            repr_latex += "{} = {}".format(lhs_str, latex_eq_rhs)
            repr_latex += r"$"

        return repr_latex



    def to_latex(self, keys=None, mode='equation'):
        """Generate list of LaTeX-formated formulas
        Upon registeration, each function should have a _repr_latex_ method.
        """
        if keys is None:
            keys = list(self.signatures.keys())

        funcs_latex = []
        for k in keys:
            # funcs_latex.append(self.func_latex(k, mode))
            func_latex = self[k]._repr_latex_()
            if mode == 'equation':
                func_latex = func_latex.strip('$')
                func_latex = r"\begin{equation}" + func_latex + r"\end{equation}"
            funcs_latex.append(func_latex)

        repr_latex = " ".join(funcs_latex)

        return beautify_latex(repr_latex).encode('utf-8').decode()

    def _repr_latex_(self):
        """Provide notebook rendering of formulas"""
        return self.to_latex()

    def detail(self):
        """Constructs a pandas dataframe from signatures"""
        return pd.DataFrame(self.signatures).T


    def simulate(self, **kwargs):
        state_funcs = []
        for name, func_key in list(self.symbol_registry.items()):
            func = self[func_key]
            update_var = getattr(func, 'update', None)
            if update_var is not None:
                state_funcs.append((func.update, func))

        return simulate(OrderedDict(state_funcs), **kwargs)

    def evaluate(self, variable, *args, **kwargs):
        """evaluates the variable

        if the variable is not present, try to parse it as a semicolon-delimited list
        """
        if not hasattr(self, variable):
            var_dict = {}
            for variable_ in variable.split(';'):
                if len(variable_.split('=')) == 2:
                    variable_name, variable_expr = variable_.strip("'").split('=')
                    var_dict[variable_name] = variable_expr
                else:
                    raise SyntaxError('cannot parse {}'.format(variable_))
            knew = from_kamodo(self, **var_dict)
            # evaluate the last variable
            result = knew.evaluate(variable=variable_name, *args, **kwargs)
            return result


        if isinstance(self[variable], np.vectorize):
            params = get_defaults(self[variable].pyfunc)
        else:
            params = get_defaults(self[variable])

        valid_arguments = valid_args(self[variable], kwargs)
        if len(params) > 0:
            if self.verbose:
                print('default parameters:', list(params.keys()))
            params.update(valid_arguments)
        else:
            if self.verbose:
                print('no default parameters')
            params = valid_arguments
            if self.verbose:
                print('user-supplied args:', list(params.keys()))
        result = eval_func(self[variable], params)
        if isinstance(result, GeneratorType):
            params = concat_solution(result, variable)
        else:
            params.update({variable: result})
        return params

    def solve(self, fprime, interval, y0,
              dense_output=True,  # generate a callable solution
              events=None,  # stop when event is triggered
              vectorized=True, ):
        from scipy.integrate import solve_ivp

        result = solve_ivp(self[fprime], interval, y0,
                           dense_output=dense_output,
                           events=events,
                           vectorized=vectorized)
        if self.verbose:
            print(result['message'])

        varname = next(iter(inspect.signature(self[fprime]).parameters.keys()))

        scope = {'result': result, 'varname': varname}

        soln_str = r"""def solution({varname} = result['t']):
            return result['sol'].__call__({varname}).T.squeeze()""".format(varname=varname)

        exec(soln_str.format(), scope)

        return scope['solution']

    def figure(self, variable, indexing='ij', return_type=False, **kwargs):
        """Generates a plotly figure for a given variable and keyword arguments"""
        result = self.evaluate(variable, **kwargs)
        signature = self.signatures[variable]
        units = signature['units']
        if units != '':
            units = '[{}]'.format(units)
        title = self.to_latex([variable], 'inline')
        title_lhs = title.split(' =')[0] + '$'
        title_short = '{}'.format(variable + units)  # something wrong with colorbar latex being vertical
        titles = dict(title=title, title_lhs=title_lhs, title_short=title_short, units=units, variable=variable)
        fig = dict()
        chart_type = None
        traces = []

        hidden_args = []
        if 'hidden_args' in self[variable].meta:
            hidden_args = self[variable].meta['hidden_args']
        arg_arrays = [result[k] for k in result if k not in hidden_args][:-1]

        arg_shapes = get_arg_shapes(*arg_arrays)

        out_dim, arg_dims = get_plot_key(result[variable].shape, *arg_shapes)

        try:
            plot_func = plot_dict[out_dim][arg_dims]['func']
        except KeyError:
            print('not supported: out_dim {}, arg_dims {}'.format(out_dim, arg_dims))
            raise

        traces, chart_type, layout = plot_func(
            result,
            titles,
            indexing=indexing,
            verbose=self.verbose,
            **kwargs)

        layout.update(
            dict(autosize=False,
                 width=700,
                 height=400,
                 margin=go.layout.Margin(
                     # l=30,
                     r=30,
                     b=32,
                     t=40,
                     pad=0),
                 ))

        fig['data'] = traces
        fig['layout'] = layout
        if return_type:
            fig['chart_type'] = chart_type
        return fig

    def plot(self, *variables, **figures):
        for k in variables:
            figures[k] = {}
        if len(figures) == 1:
            variable, kwargs = list(figures.items())[0]
            fig = self.figure(variable, return_type=True, **kwargs)
            if fig['chart_type'] is None:
                raise AttributeError("No chart_type for this trace")
            else:
                if self.verbose:
                    print('chart type:', fig['chart_type'])
                return go.Figure(data=fig['data'], layout=fig['layout'])
        else:
            traces = []
            layouts = []
            for variable, kwargs in list(figures.items()):
                fig = self.figure(variable, **kwargs)
                traces.extend(fig['data'])
                layouts.append(fig['layout'])
            # Todo: merge the layouts instead of selecting the last one
            return go.Figure(data=traces, layout=layouts[-1])


class KamodoAPI(Kamodo):
    """JSON api wrapper for kamodo services"""
    def __init__(self, url_path, **kwargs):
        self._url_path = url_path
        super(KamodoAPI, self).__init__(**kwargs)

        self._kdata = self._get(self._url_path)

        self._defaults = {}
        self._data = {}

        for k, v in self._kdata.items():
            self[v['lhs']] = kamodofy(units=v['units'])

            # get defaults for this func
            default_path = '{}/{}/defaults'.format(self._url_path, k)
            self._defaults[k] = self._get(default_path)

            # get cached data (result of calling with no args)
            data_path = '{}/{}/data'.format(self._url_path, k)
            self._data[k] = self._get(data_path)

            # get kamodofied function
            self[k] = kamodofy(
                self.load_func(k),
                data=self._data[k],
                units=v['units'])

    def _get(self, url_path):
        req_result = requests.get(url_path)
        try:
            result = req_result.json() # returns a dictionary maybe
        except json.JSONDecodeError as m:
            print('could not decode request {}'.format(req_result.text))
            raise json.JSONDecodeError(m)

        if isinstance(result, str):
            result_dict = json.loads(result, object_hook=deserialize)
        else:
            result_dict = result
        return result_dict

    def _call_func(self, func_name, **kwargs):
        """construct the url and call api with params"""
        url_path = '{}/{}'.format(self._url_path, func_name)
        if self.verbose:
            print('querying {}'.format(url_path))
        params = []
        for k, v in kwargs.items():
            params.append((k, json.dumps(v, default=serialize)))
        result = requests.get(
            url=url_path,
            params=params).json() #returns a dictionary

        if isinstance(result, str):
            if self.verbose:
                print('reloading json as str')
            result = json.loads(result, object_hook=deserialize)
        else:
            result = deserialize(result)
        return result

    def load_func(self, func_name):
        """loads a function signature"""
        signature = []
        for arg, arg_default in self._defaults[func_name].items():
            signature.append(forge.arg(arg, default=arg_default))

        @forge.sign(*signature)
        def api_func(*args, **kwargs):
            """API function"""
            return self._call_func(func_name, **kwargs)

        api_func.__name__ = func_name
        api_func.__doc__ = "{}/{}".format(self._url_path, func_name)
        return api_func

def compose(**kamodos):
    """Kamposes multiple kamodo instances into one"""
    kamodo = Kamodo()
    for kname, k in kamodos.items():
        for name, symbol in k.symbol_registry.items():
            signature = k.signatures[name]
            meta = k[symbol].meta
            data = getattr(k[symbol], 'data', None)

            rhs = signature['rhs']
            registry_name = '{}_{}'.format(name, kname)
            if (rhs is None) | hasattr(rhs, '__call__'):
                kamodo[registry_name] = kamodofy(k[symbol], data=data, **meta)
            else:
                kamodo[registry_name] = str(rhs)

    return kamodo

def from_kamodo(kobj, **funcs):
    """copies a kamodo object, inserting additional functions"""
    knew = Kamodo()
    for name, signature in kobj.signatures.items():
        symbol = signature['symbol']
        knew[symbol] = kobj[symbol]
    for symbol, func in funcs.items():
        knew[symbol] = func
    return knew
