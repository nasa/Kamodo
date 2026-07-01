# -*- coding: utf-8 -*-
"""
Copyright © 2017 United States Government as represented by the Administrator, National Aeronautics and Space Administration.
No Copyright is claimed in the United States under Title 17, U.S. Code.  All Other Rights Reserved.
"""

import asyncio
import os
import socket
import copy
import time
import functools
import inspect
import itertools
import json
import re

import types
from collections import OrderedDict
from collections import UserDict
from inspect import getfullargspec
from types import GeneratorType

import forge
import pandas as pd

import yaml
import numpy as np
import plotly.graph_objs as go
import requests
from plotly.subplots import make_subplots

import sympy
from sympy import Eq, Pow
from sympy import Expr
from sympy import Symbol, symbols, Function
from sympy import lambdify
from sympy import latex
from sympy.abc import _clash
from sympy.core.function import UndefinedFunction
from sympy.parsing.latex import parse_latex
from sympy.physics import units as sympy_units
from sympy.physics.units import Dimension

from plotting import get_ranges
from plotting import plot_dict, get_arg_shapes, symbolic_shape
from rpc.proto import Server, ssl
from rpc.proto import capnp, KamodoRPC, FunctionRPC, kamodo_capnp, rpc_map_to_dict
from rpc.proto import from_rpc_literal, to_rpc_literal, to_rpc_expr
from rpc.proto import wrap_async
# from util import to_arrays, cast_0_dim
from util import beautify_latex
from util import concat_solution
from util import construct_signature

from util import get_arg_units
from util import get_defaults, valid_args, eval_func
from util import get_dimensions
from util import kamodofy
from util import np
from util import partial
from util import reserved_names
from util import serialize, deserialize
from util import sign_defaults
from util import simulate
from util import unify, get_abbrev, get_expr_unit
from util import unit_subs

# try:
#     from sympy.parsing.sympy_parser import parse_expr
# except ImportError:  # may occur for certain versions of sympy
#     from sympy import sympify as parse_expr
# from IPython.display import Latex
# -

_clash['rad'] = Symbol('rad')
_clash['deg'] = Symbol('deg')


def parse_expr(*args, **kwargs):
    try:
        return sympy.sympify(*args, **kwargs)
    except:
        raise NameError('cannot parse {}, {}'.format(args, kwargs))


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
    unit_expr = parse_expr(unit_str.replace('^', '**'), locals=_clash)
    try:
        unit = unit_expr.subs(units)
    except:
        raise NameError(
            'something wrong with unit str [{}], type {}, {}'.format(
                unit_str, type(unit_str), unit_expr))
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



def alphabetize(symbols):
    """alphabetical ordering for the input symbols"""
    # can't use > on symbols, so need to convert to str first
    return tuple(
        Symbol(symbol_) for symbol_ in sorted([str(_) for _ in symbols]))


def reorder_symbol(defaults, default_non_default_parameter, symbol):
    """
    Changed the value of symbol based on the order of
    default_non_default_parameter, if the order is not appropriate then
    reorder both default and non default  parameter alphabetically,
    updates symbol value and make sure default comes always after non default
    parameter
    """
    default_parameter_keys = set(defaults.keys())
    default_non_default_parameter = set(default_non_default_parameter)
    non_default_parameters = default_non_default_parameter - default_parameter_keys
    non_default_parameters = sorted(non_default_parameters)
    non_default_parameters_prev = non_default_parameters
    default_parameters = sorted(default_parameter_keys)
    non_default_parameters.extend(default_parameters)
    new_symbol = tuple(Symbol(symbol_) for symbol_ in
                       [str(_) for _ in non_default_parameters])
    dict_symbol = {}
    for i in range(len(new_symbol)):
        dict_symbol["var{0}".format(i)] = new_symbol[i]

    new_symbol = Function(symbol.name)(*dict_symbol.values())
    if new_symbol.args != symbol.args:
        symbol = new_symbol
    return symbol


def default_inheritance_check(rhs_expr, lhs_expr):
    """
    Checking  if default parameter comes first and then non-default
    parameter, ex:G(X,Y).In such cases length of rhs_expr will be greater
    than 1, 1st and 2nd argument of rhs_expr will be sympy Symbol and
    sympy function data type.If condition satisfies, then checking if the
    first argument for both rhs_expr and lhs_expr is same or not.If not same
    then raise syntax error.Same check will happen in case more than 1 default
    value, additionally it will discard few cases to pass all the other
    functionalities
    """
    try:
        if isinstance(rhs_expr.args[0], Symbol):
            if len(rhs_expr.args) == 2:
                if type(rhs_expr.args[1]) != rhs_expr.args[1].name:
                    if rhs_expr.args[0] != lhs_expr.args[0] and not \
                            isinstance(rhs_expr.args[1], Symbol):
                        raise SyntaxError('Ordering error')
            elif len(rhs_expr.args) > 2 and rhs_expr.args != lhs_expr.args \
                    and not len(rhs_expr.args[1].args) > 0:
                for each in rhs_expr.args:
                    if isinstance(each, Symbol) or isinstance(each, Pow):
                        pass
                    else:
                        if rhs_expr.args[0] != lhs_expr.args[-1]:
                            raise SyntaxError(f'Ordering error {lhs_expr} = {rhs_expr}')
    except IndexError:
        pass
    except TypeError:
        pass
    except AttributeError:
        pass


def dimensionless_unit_check(sym_name, arg_units):
    """
    Check if args_units has dimensionless item using regex,
    if it has then 'ordered_unit' will return blank item  for the
    corresponding key and will update 'arg_units' with blank-""
    """
    try:
        sym_name_split = sym_name.split(')')
        sym_name_split.pop(-1)
        for i in sym_name_split:
            if '(' in i and '[' in i:
                flag = True
            else:
                flag = False
        if flag:
            lhs_args_temp = re.findall(r"\((.+)\)", sym_name)[0]
            ordered_unit = re.findall(
                r"((?<!\])|(?<=\[)[^\[\],]*)\]?(?:,|\)|$)",
                lhs_args_temp)
            if "" in ordered_unit:
                counter = 0
                for k, v in arg_units.items():
                    arg_units[k] = ordered_unit[counter]
                    counter = counter + 1
        return arg_units

    except KeyError:
        pass


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
        return False  # no paren hence no arguments
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
            extract_units('f')
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

        try:
            if len(args) != len(arg_units) and len(arg_units) > 0:
                for i in range(len(args)):
                    if '[' not in args[i]:
                        arg_units.insert(i, "")
        except IndexError:
            pass
        except NotImplementedError:
            pass

        for arg, unit in zip(args, arg_units):
            unit_dict[arg.replace('[{}]'.format(unit), '')] = unit
        if len(all_units) == len(arg_units) and "" not in arg_units:
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
    return symbols(
        [a for a in getfullargspec(func).args if a not in hidden_args])


class Kamodo(UserDict):
    """Kamodo base class demonstrating common API for scientific resources.

    This API provides access to scientific fields and their properties through:
    
    * Interpolation of variables at user-defined points
    * Automatic unit conversions
    * Function composition convenient for coordinate transformations and data pipelining

    Note: While intended for with space weather applications, the kamodo base class
    was designed to be as generic as possible, and should be applicable to a wide range
    of scientific domains and disciplines.
    """

    def __init__(self, *funcs, **kwargs):
        """Initialize Kamodo object with functions or by keyword
        
        ** Inputs **

        * ** funcs ** - *(optional)* list of (str) expressions to register in f(x)=x format

        * ** kwargs ** - *(optional)* key,value pairs of functions to register
            * key - left-hand-side symbol (str)
            * value - can be one of:
                * latex or python (str) expression e.g. "x^2-x-1"
                * kamodofied function with appropriate .meta and .data attributers (see [@kamodofy](#kamodofy) decorator)
                * lambda function (having no meta or data attributes)

        * ** verbose ** - *(optional)* (`default=False`) flag to turn on all debugging print statements


        ** returns ** - dictionary-like kamodo object of (symbol, function) pairs

        usage:

        ```python
            kobj = Kamodo(
                'f(x[cm])[kg/m^3]=x^2-x-1', # full expressions with units
                area = kamodofy(lambda x: x*x, units='cm^2'), # kamodofied functions
                h = 'sin(x)', # key-value expressions
                )
        ```

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
            print('symbol arguments: >>> {} {} <<<'.format(symbol.args,
                                                           type(symbol.args)))
            print('free_symbols: >>> {} {} <<<'.format(free_symbols,
                                                       type(free_symbols)))
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
                    raise NameError(
                        "Mismatched symbols {} and {}".format(lhs_args,
                                                              free_symbols))
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
                print(
                    'lambda {} = {} labmdified with numexpr'.format(symbol.args,
                                                                    rhs_expr))
        except:  # numexpr not installed
            func = lambdify(symbol.args, rhs_expr,
                            modules=['numpy', composition])

            if self.verbose:
                print(
                    'lambda {} = {} lambdified with numpy and composition:'.format(
                        symbol.args, rhs_expr))
                for k, v in composition.items():
                    print('\t', k, v)
        signature, defaults = sign_defaults(symbol, rhs_expr, composition)

        return signature(func)

    def update_unit_registry(self, func, arg_units):
        """Inserts unit functions into registry"""
        lhs, unit_dict = extract_units(func)
        if self.verbose:
            print('extracted lhs units: {}'.format(lhs))
            for k, v in unit_dict.items():
                print('  ', k, v)
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

    def register_signature(self, symbol, units, lhs_expr, rhs_expr, arg_units):
        # if isinstance(units, str):
        unit_str = units
        if self.verbose:
            print('unit str {}'.format(unit_str))

        if rhs_expr is None:
            lambda_ = symbols('lambda', cls=UndefinedFunction)
            rhs_expr = lambda_(*symbol.args)

        self.signatures[str(type(symbol))] = dict(
            symbol=symbol,
            units=unit_str,
            lhs=lhs_expr,
            rhs=rhs_expr,
            arg_units=arg_units
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

        lhs_symbol = self.check_or_replace_symbol(lhs_symbol, rhs_args)
        units = lhs_units
        if hasattr(func, 'meta'):
            if self.verbose:
                print('function has metadata', func.meta.keys())
            if len(lhs_units) > 0:
                if lhs_units != func.meta['units']:
                    raise NameError('Units mismatch:{} != {}'.format(lhs_units,
                                                                     func.meta[
                                                                         'units']))
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
            unit_dict = {arg: get_unit(arg_unit) for arg, arg_unit in
                         arg_units.items()}
            unit_expr = lhs_symbol.subs(unit_dict)
            self.unit_registry[lhs_symbol] = unit_expr
            self.unit_registry[unit_expr] = get_unit(units)
        else:
            self.unit_registry[lhs_symbol] = get_unit(units)
        self.register_signature(lhs_symbol, units, lhs_expr, rhs, arg_units)
        try:
            if hasattr(func, '_repr_latex_'):
                symbol = list(self.signatures.items())[0][1]['symbol']
                func._repr_latex_ = lambda: self.func_latex(
                    str(type(lhs_symbol)),
                    mode='inline')
                key = list(self.signatures.keys())[0]
                self.signatures[key]['rhs'] = func._rhs_
            else:
                func._repr_latex_ = lambda: self.func_latex(str(type(
                    lhs_symbol)), mode='inline')

        except AttributeError:
            # happens with bound methods
            pass
        super(Kamodo, self).__setitem__(lhs_symbol, func)  # assign key 'f(x)'
        super(Kamodo, self).__setitem__(type(lhs_symbol),
                                        self[lhs_symbol])  # assign key 'f'
        self.register_symbol(lhs_symbol)

    def __setitem__(self, sym_name, input_expr):
        """Assigns a function or expression to a new symbol, performing
        automatic function composition and inserting unit conversions where appropriate.

        * ** sym_name ** - function symbol to associate with right-hand-side in one of the following formats:
            - f - a lone fuction symbol (alphabetic argument ordering)
            - f(z,x,y) - explicit argument ordering
            - f[kg] - output unit assignment
            - f(x[cm])[kg] - output and input unit assignment

        * ** input_expr ** - rhs string or kamodofied function, one of:
            * right-hand-side expression: python or latex str (e.g.`x^2-x-1`)
            * kamodofied function with appropriate .meta and .data attributers (see [@kamodofy](#kamodofy))
            * lambda function (having no meta or data attributes)

        Raises:
            - NameError when left-hand-side units incompatible with right-hand-side expression
        
        returns: None

        usage:

        Setting left-hand-side units will automatically trigger unit conversion
        
        ```py
        kobj = Kamodo()
        kobj['radius[m]'] = 'r'
        kobj['area[cm^2]'] = 'pi * radius^2'
        kobj
        ```

        The above `kobj` will render in a Jupyter notebook like this:

        $$\\operatorname{radius}{\\left(r \\right)}[m] = r$$
        
        $$\\operatorname{area}{\\left(r \\right)}[cm^{2}] = 10000 \\pi \\operatorname{radius}^{2}{\\left(r \\right)}$$

        Kamodo will raise an error if left-hand-side units are incompatible with the right-hand-side expression
        
        ```py
        kobj = Kamodo()
        kobj['area[cm^2]'] = 'x^2' # area has units of cm^2
        try:
            kobj['g(x)[kg]'] = 'area' # mass not compatible with square length
        except NameError as m:
            print(m)
        ```
        
        output:
        
        $$\\text{cannot convert area(x) [centimeter**2] length**2 to g(x)[kilogram] mass}$$
        

        """
        if not isinstance(sym_name, str):
            sym_name = str(sym_name)
        symbol, args, lhs_units, lhs_expr = self.parse_key(sym_name)
        if hasattr(input_expr, '__call__'):
            input_expr = copy_func(input_expr)
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

            default_inheritance_check(rhs_expr, lhs_expr)
            if not isinstance(symbol, Symbol):
                if isinstance(lhs_expr, Symbol):
                    symbol = Function(lhs_expr)(*tuple(alphabetize(
                        rhs_expr.free_symbols)))
                else:  # lhs is already a function
                    symbol = lhs_expr
                lhs_str = str(symbol)
                sym_name = sym_name.replace(str(lhs_expr), lhs_str)
            if self.verbose:
                print('unit registry contents:')
                for k, v in self.unit_registry.items():
                    print('\t', k, type(k), v)
            if '[' in sym_name:
                if self.verbose:
                    print(
                        'updating unit registry with {} -> {}'.format(sym_name,
                                                                      rhs_expr))
                rhs = rhs_expr
                arg_units = get_arg_units(rhs_expr, self.unit_registry)
                if self.verbose:
                    print(arg_units)

                sym_name_bkup = sym_name

                sym_name = self.update_unit_registry(sym_name, arg_units)
                if self.verbose:
                    print('unit registry update returned', sym_name,
                          self.unit_registry.get(symbol))
            else:

                if self.verbose:
                    print(sym_name,
                          symbol,
                          'had no units. Getting units from {}'.format(
                              rhs_expr))
                expr_unit = get_expr_unit(rhs_expr, self.unit_registry,
                                          self.verbose)
                arg_units = get_arg_units(rhs_expr, self.unit_registry)

                if self.verbose:
                    print('registering {} with {} {}'.format(symbol, expr_unit,
                                                             arg_units))

                if (symbol not in self.unit_registry) and (
                        expr_unit is not None):
                    self.unit_registry[symbol] = symbol.subs(arg_units)
                    self.unit_registry[symbol.subs(arg_units)] = expr_unit

                if expr_unit is not None:
                    expr_dimensions = Dimension(get_dimensions(expr_unit))
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
                print('unit registry to resolve units:')
                for k, v in self.unit_registry.items():
                    print('\t{}:{}'.format(k, v))

            units = get_expr_unit(symbol, self.unit_registry)
            if Dimension(get_dimensions(units)) != Dimension(1):
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

            for k, v in arg_units.items():
                if str(v) == 'Dimension(1)':
                    arg_units[k] = ''

            signature, defaults = sign_defaults(symbol, rhs_expr, composition)
            default_non_default_parameter = []
            try:
                for parm in signature.parameters:
                    default_non_default_parameter.append(parm.name)
            except KeyError:
                pass

            # symbol = reorder_symbol(defaults, default_non_default_parameter,
            #                         symbol)

            if len(defaults) > 0:
                symbol = reorder_symbol(defaults, default_non_default_parameter,
                                        symbol)

            try:
                arg_units = dimensionless_unit_check(sym_name_bkup,
                                                     arg_units)
            except UnboundLocalError:
                if len(rhs.args) > 0:
                    try:
                        split_rhs_args = str(rhs.args).split(',')[0]
                        if ('(' in split_rhs_args) and (')' not in
                                                        split_rhs_args):
                            arg_units = {}
                            units = ""
                    except IndexError:
                        pass

            meta = dict(units=units, arg_units=arg_units)
            func.meta = meta
            func.data = None
            self.register_signature(symbol, units, lhs_expr, rhs_expr,
                                    arg_units)
            func._repr_latex_ = lambda: self.func_latex(str(type(symbol)),
                                                        mode='inline')

            self.register_symbol(symbol)
            func._preserved_ = lambda: self.func_latex(str(type(symbol)),
                                                       mode='inline')
            func._rhs_ = list(self.signatures.items())[0][1]['rhs']

            super(Kamodo, self).__setitem__(symbol, func)
            super(Kamodo, self).__setitem__(type(symbol), self[symbol])

    def __getitem__(self, key):
        """Given a symbol string, retrieves the corresponding function.

        input: **key** - string or function symbol

    
        ** returns**: the associated function

        ** usage **:

        Rretrieval by function name:

        ```python
            kobj['f'] = 'x^2-x-1'
            f = kobj['f']
            f(3) # returns 5
        ```

        It is also possible to retreive by function symbol:

        ```python
            from kamodo import sympify
            fsymbol = sympify('f') # converts str to symbol

            kobj['f'] = 'x^2-x-1'
            f = kobj[fsymbol]
            f(3) # returns 5
        ```

        """
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
        """

        Retrieves a given function as an attribute.
        
        **input** - **name** of function to retrieve

        **returns** the associated function

        Usage:

        ```py
        k = Kamodo(f='x^2-x-1')
        k.f
        ```
        The above renders as follows in a jupyter notebook

        $f{\\left(x \\right)} = x^{2} - x - 1$

        """
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    def __delattr__(self, name):

        if name in self:
            symbol = self.symbol_registry[name]
            if self.verbose:
                print('__delattr__: removing name {}, symbol {}'.format(name,
                                                                        symbol))
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
            if self.verbose:
                print(
                    '__delitem__: removing {} from symbol_registry'.format(key))
            symbol = self.symbol_registry.pop(key, None)
            if symbol is None:
                symbol, args, unit_dict, parsed = parse_lhs(
                    key, self.symbol_registry, self.verbose)
        else:
            symbol = key
            if self.verbose:
                print('__delitem__: received non str key {}'.format(symbol))
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
            self.signatures.pop(str(key_), None)
            self.signatures.pop(str(type(key_)), None)
            self.unit_registry.pop(key_, None)


    def func_latex(self, key, mode='equation'):
        """get a latex string for a given function key"""
        repr_latex = ""
        lhs = self.signatures[key]['symbol']
        rhs = self.signatures[key]['rhs']
        units = self.signatures[key]['units']
        arg_units = get_arg_units(lhs, self.unit_registry)

        for k, v in arg_units.items():
            if str(v) == 'Dimension(1)':
                arg_units[k] = ''

        if len(units) > 0:
            units = '{}'.format(get_abbrev(units))
        else:
            units = ''

        if len(arg_units) > 0:
            arg_strs = []
            for arg, arg_unit in arg_units.items():
                arg_strs.append("{}[{}]".format(
                    latex(arg),
                    latex(get_abbrev(arg_unit),
                          fold_frac_powers=True,
                          fold_short_frac=True,
                          root_notation=False,
                          )))
            lhs_str = "{}({})".format(latex(type(lhs)), ",".join(arg_strs))
        else:
            lhs_str = latex(lhs)
            # lhs_str = "{}({})".format(
            #     latex(type(lhs)),
            #     ','.join([latex(s) for s in lhs.args]))
        if len(units) > 0:

            lhs_str += "[{}]".format(
                latex(parse_expr(units.replace('^', '**'), locals=_clash),
                      fold_frac_powers=True,
                      fold_short_frac=True,
                      root_notation=False,
                      ))
            dimension_less_args = lhs_str.find('[]')
            if dimension_less_args != -1:
                lhs_str = lhs_str.replace('[]', '')

        latex_eq = ''
        latex_eq_rhs = ''
        if isinstance(rhs, str):
            latex_eq_rhs = rhs
        elif hasattr(rhs, '__call__') | (rhs is None):
            lambda_ = symbols('lambda', cls=UndefinedFunction)
            # latex_eq = latex(Eq(lhs, lambda_(*lhs.args)), mode=mode)
            latex_eq_rhs = latex(lambda_(*lhs.args))  # no $$
        else:
            # latex_eq = latex(Eq(lhs, rhs), mode=mode)
            latex_eq_rhs = latex(rhs)  # no $$ delimiter

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

        ** inputs **:

        * keys - (optional) list(str) of registered functions to generate LaTeX from

        * mode - (optional) string determines to wrap formulas
            * 'equation' (default) wraps formulas in `begin{equation} ... end{equation}`

            * 'inline': wraps formulas in dollar signs

        ** returns **: LaTeX-formated string
        
        Note: This function does not need to be called directly for rendering in jupyter
        because the _repr_latex_ method is automatically attached.

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
        """Provides notebook rendering of kamodo object's registered functions.

        ** inputs ** - N/A

        ** returns ** latex string - obtained  from `to_latex` method

        Usage:

        ```python
        k = Kamodo(f='x^2-x-1')
        k
        ```

        When placed on a line by itself, the above object will be rendered by jupyter notebooks like this:

        \\begin{equation}f{\\left(x \\right)} = x^{2} - x - 1\\end{equation}

        More on the _repr_latex_ method can be found [here](https://ipython.readthedocs.io/en/stable/api/generated/IPython.display.html) 
        
        Note: Registered functions are also equiped with their own `_repr_latex_` method.

        """
        return self.to_latex()

    def detail(self):
        """Constructs a pandas dataframe from signatures

        ** inputs ** - N/A

        ** returns ** - pandas dataframe

        usage:

        ```python
        k = Kamodo('rho(x[cm])[g/cm^3]=x^2', g = 'x+y')
        k.detail()
        ```
        outputs:


        <table border="1" class="dataframe">  <thead>    <tr style="text-align: right;">      <th></th>      <th>symbol</th>      <th>units</th>      <th>lhs</th>      <th>rhs</th>      <th>arg_units</th>    </tr>  </thead>  <tbody>    <tr>      <th>rho</th>      <td>rho(x)</td>      <td>g/cm**3</td>      <td>rho(x)</td>      <td>x**2</td>      <td>{\'x\': \'cm\'}</td>    </tr>    <tr>      <th>g</th>      <td>g(x, y)</td>      <td></td>      <td>g</td>      <td>x + y</td>      <td>{}</td>    </tr>  </tbody></table>
        """


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
        """Evaluate a given function variable using kwargs.

        If the variable is not present, try to parse it as an equation and evaluate the expression.

        ** inputs **:

        * variable - str:
            * function string name to evaluate
            * semicolon delmitted list of equations, the last of which will be evaluated
        * args - not presently used
        * kwargs - key-word arguments passed to function (required)

        ** returns **: dictionary of input kwargs and output {variable: self.variable(**kwargs)}

        ** usage **:
        
        ```py
        k = Kamodo(f='x+y')

        result = k.evaluate('f', x=3, y=4)['f']
        assert k.f(3,4) == result
        assert k.evaluate('g=f+3', x=3, y=4)['g'] == result+3
        assert k.evaluate('g=f+3;h=g+2', x=3, y=4)['h'] == result+3+2
        ```

        """
        if not hasattr(self, variable):
            var_dict = {}
            for variable_ in variable.split(';'):
                if len(variable_.split('=')) == 2:
                    variable_name, variable_expr = variable_.strip("'").split(
                        '=')
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
            if (len(params) == 1) & isinstance(list(params.values())[0],
                                               GeneratorType):
                if self.verbose:
                    print('found parametric generator')
                # raise NotImplementedError('Animations not yet supported')
                params.update({variable: result})
            else:
                if self.verbose:
                    print(
                        'evaluate found generator function. Concatonating solutions')
                    print('generator params:', params)

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
            return result['sol'].__call__({varname}).T.squeeze()""".format(
            varname=varname)

        exec(soln_str.format(), scope)

        return scope['solution']

    def to_rpc_meta(self, key):
        """create rpc metadata"""
        meta = self[key].meta

        units = meta.get('units')
        if units is None:
            units = ''

        arg_unit_entries = []
        arg_units = meta.get('arg_units')  # may be None
        if arg_units is not None:
            for k, v in arg_units.items():
                arg_unit_entries.append({'key': k, 'value': v})

        citation = meta.get('citation')
        if citation is None:
            citation = ''

        equation = meta.get('equation')
        if equation is None:
            equation = latex(self.signatures[key]['rhs'])

        if self.verbose:
            print('equation: {}'.format(equation))

        hidden_args = meta.get('hidden_args')
        if hidden_args is None:
            hidden_args = []

        return kamodo_capnp.Kamodo.Meta(
            units=units,
            argUnits=dict(entries=arg_unit_entries),
            citation=citation,
            equation=equation,
            hiddenArgs=hidden_args,
        )

    def register_rpc_field(self, key):
        func = self[key]
        signature = self.signatures[key]
        field = kamodo_capnp.Kamodo.Field.new_message(
            func=FunctionRPC(func),
            meta=self.to_rpc_meta(key),
        )
        self._kamodo_rpc[key] = field

    def serve(self, host='localhost', port='60000', certfile=None, keyfile=None):
        """Serve registered functions using Kamodo-RPC spec

        Uses asyncio and capnp proto 

        ** inputs **:

        * host - str: localhost, ipv4 or ipv6 address (localhost by default)
        * port - str: port to serve from
        * certfile - certficicate to authenticate clients
        * keyfile - private key file to authenticate

        ** returns **: None

        ** usage **:
        
        ```py
        k = Kamodo(f='x+y')

        k.serve() # start rpc server
        ```

        see kamodo/rpc/kamodo.capnp
        """
        self._kamodo_rpc = KamodoRPC()
        for key in self.signatures:
            if self.verbose:
                print('serving {}'.format(key))
            self.register_rpc_field(key)
        if self.verbose:
            print(f'serving with \n {certfile}\n {keyfile}')
        self.async_server = Server(self._kamodo_rpc)
        asyncio.run(self.async_server.serve(host, port, certfile, keyfile))

    def figure(self, variable, indexing='ij', **kwargs):
        """Generates a plotly figure for a single variable and keyword arguments
        
        ** inputs **:

        * variable: the name of a previously registered function
        * kwargs: {arg: values} to pass to registered function
        * indexing: determines order by which 2d matrices are given (affects contour_plot, carpet_plot, and plane)

        ** returns **: plotly [figure](https://plotly.com/python/figure-structure/) (dict-like)

        raises: SyntaxError if variable not found
        """
        result = self.evaluate(variable, **kwargs)
        signature = self.signatures[variable]
        units = signature['units']
        if units != '':
            units = '[{}]'.format(units)
        title = self.to_latex([variable], 'inline')
        title_lhs = title.split(' =')[0] + '$'
        title_short = '{}'.format(
            variable + units)  # something wrong with colorbar latex being vertical
        titles = dict(
            title=title, title_lhs=title_lhs,
            title_short=title_short, units=units, variable=variable)
        fig = dict()
        chart_type = None
        traces = []

        hidden_args = []
        if 'hidden_args' in self[variable].meta:
            hidden_args = self[variable].meta['hidden_args']
        arg_arrays = [result[k] for k in result if k not in hidden_args][:-1]

        arg_shapes = get_arg_shapes(*arg_arrays)

        if isinstance(result[variable], GeneratorType):
            # if evaluate returned a generator, then assume animation
            raise NotImplementedError('Animations not yet supported!')

        try:
            out_dim, *arg_dims = symbolic_shape(result[variable].shape,
                                                *arg_shapes)
        except IndexError:
            print('could not interpret shapes from variable shape {}'.format(
                result[variable].shape))
            print('argument shapes: ', *arg_shapes)
            raise

        try:
            plot_func = plot_dict[out_dim][tuple(arg_dims)]['func']
        except KeyError:
            print('not supported: out_dim {}, arg_dims {}'.format(out_dim,
                                                                  arg_dims))
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
        try:
            if len(signature['arg_units']) == 1:
                try:
                    args_unit = signature['arg_units']
                    x_axis_unit = \
                        [args_unit[i][0] for i in sorted(args_unit.keys())][0]
                    y_axis_unit = signature['units']
                    x_last_index = layout.xaxis.title.text.rindex('$')
                    y_last_index = layout.yaxis.title.text.rindex('$')
                    new_xaxis_label = layout.xaxis.title.text[
                                      :x_last_index] + ' ' + f'[{x_axis_unit}]' + \
                                      layout.xaxis.title.text[x_last_index:]
                    new_yaxis_label = layout.yaxis.title.text[
                                          y_last_index] + str(
                        signature['symbol'].name) + ' ' + f'[{y_axis_unit}]' + \
                                      layout.yaxis.title.text[y_last_index]

                    layout['xaxis']['title']['text'] = new_xaxis_label
                    layout['yaxis']['title']['text'] = new_yaxis_label
                except AttributeError:
                    pass
                except IndexError:
                    pass
            elif len(signature['arg_units']) == 2:
                try:
                    arg_units = signature['arg_units'].values()
                    arg_keys = list(signature['arg_units'].keys())
                    x_ax_unit = f" [{list(arg_units)[0]}]"
                    y_ax_unit = f" [{list(arg_units)[1]}]"
                    dolr_char = layout.xaxis.title.text.rindex('$')
                    new_xaxis_label = layout.xaxis.title.text[dolr_char] + ' ' + \
                                      arg_keys[0] + " " + f"{x_ax_unit}" + \
                                      layout.xaxis.title.text[dolr_char]
                    new_yaxis_label = layout.xaxis.title.text[dolr_char] + ' ' + \
                                      arg_keys[1] + " " + f"{y_ax_unit}" + \
                                      layout.xaxis.title.text[dolr_char]

                    layout['xaxis']['title']['text'] = new_xaxis_label
                    layout['yaxis']['title']['text'] = new_yaxis_label

                except AttributeError:
                    pass
                except IndexError:
                    pass
        except TypeError:
            pass

        fig['data'] = traces
        fig['layout'] = layout
        return go.Figure(fig).update_traces(meta=chart_type)
        # if return_type:
        #     fig['chart_type'] = chart_type
        # return fig

    def plot(self, *variables, plot_partial={}, **figures):
        """Generates a plotly figure from multiple variables and keyword arguments

        ** inputs **:

        * variable: the name of a previously registered function
        * plot_partial: dict(dict) of {varname: {arg: values}} partial arguments to fix
        * figures: dict {variable: {arg: values}} to pass to registered function
        
        ** returns **: plotly [figure](https://plotly.com/python/figure-structure/) (dict-like).
        When run in a jupyter notebook, an inline plotly figure will be displayed.

        ** raises **:

        * TypeError when required function arguments are not specified
        * KeyError when no plotting function can be associated with input/output shapes
        
        ** usage **:
        
        ```python
        k = Kamodo(
            f=lambda x=np.array([2,3,4]): x**2-x-1,
            g='sin(x)')

        k.plot('f') # plots f using default arguments for f

        k.plot(f={x:[3,4,5]}, g={x{-2, 3, 4}}) # plots f at x=[3,4,5] and g at [-2,3,4]
        ```

        Use the `plot_partial` keyword to lower the dimensionality of a function by fixing some of its variables:

        ```python
        from kamodo import kamodofy, gridify, Kamodo
        from scipy.interpolate import RegularGridInterpolator
        import numpy as np

        # define sample coordinate and data arrays
        t, lon, lat, ht = np.linspace(0.,24.,10),
            np.linspace(0.,360.,20),
            np.linspace(-90.,90.,50),
            np.linspace(100.,10000.,250)

        variable = np.reshape(
            np.linspace(0., np.pi, 2500000), # data
            (10, 20, 50, 250)) # match coordinate arrays

        # define and kamodofy interpolating function
        rgi = RegularGridInterpolator((t, lon, lat, ht), variable, bounds_error = False, fill_value=np.NaN)

        @kamodofy(units='m/s', data=variable)
        def interpolator(xvec):
            return rgi(xvec)

        #gridify same function
        interpolator_grid = kamodofy(gridify(interpolator, time = t, lon=lon, lat = lat,  height = ht), units='m/s', data=variable, arg_units={'time':'hr','lon':'deg','lat':'deg','height':'km'})

        #register in a new kamodo object
        kamodo_object = Kamodo()
        kamodo_object['v'] = interpolator
        kamodo_object['v_ijkl'] = interpolator_grid # 4 dimensional
        kamodo_object

        kamodo_object.plot(v_ijkl = {'time':10.,'lat':90.})
        ```

        The above line raises an Error: not supported: out_dim ('N', 'M'), arg_dims [(1,), ('N',), (1,), ('M',)]

        Since there is no straight-foward way to plot a high-dimensional function, we can use `plot_partial` instead
        to get a lower-dimensional slice instead: 

        ```python
        kamodo_object.plot('v_ijkl', plot_partial={'v_ijkl': {'lon':10.,'lat':90.}})
        ```

        See also [@partial](#partial) decorator to fix a function's arguments at registration time.

        """
        if len(plot_partial) > 0:
            # kpartial = from_kamodo(self)  # copy kamodo object
            # kpartial = {}
            kpartial = Kamodo()
            for k, v in plot_partial.items():
                regname = k
                kpartial[regname] = partial(self[k], **v)
            return kpartial.plot(*variables, **figures)

        for k in variables:
            figures[k] = {}
        if len(figures) == 1:
            variable, kwargs = list(figures.items())[0]
            fig = self.figure(variable, **kwargs)
            # if fig['chart_type'] is None:
            #     raise AttributeError("No chart_type for this trace")
            # else:
            #     if self.verbose:
            #         print('chart type:', fig['chart_type'])
            return go.Figure(
                data=fig['data'],
                layout=fig['layout'])
        else:
            traces = []
            layouts = []
            for variable, kwargs in list(figures.items()):
                fig = self.figure(variable, **kwargs)
                traces.extend(fig['data'])
                layouts.append(fig['layout'])

            return go.Figure(data=traces, layout=layouts[-1])


class KamodoClient(Kamodo):
    def __init__(self, host='localhost', port='60000', certfile="selfsigned.cert", **kwargs):
        """CapnProto Kamodo client
        
        Abstracts a remote kamodo server using capn proto binary message types

        ** inputs **:

        * host - str: localhost, ipv4 or ipv6 address (localhost by default)
        * port - str: port to serve from
        * certfile - certficicate to authenticate clients (selfsigned.cert)

        ** returns **: Kamodo object with server-side functions

        ** usage **:
        
        ```py
        k = KamodoClient() # connect to localhost:60000 by default
        k.f # assuming f is registered on remote
        ```

        """
        super(KamodoClient, self).__init__(**kwargs)
        self._expressions = {}  # expressions for server-side pipelining
        self._rpc_funcs = {}
        self.host = host  # rpc functions (may be served to downstream applications)
        self.port = port
        self.certfile = certfile
        if host and port is not None:
            self.connect(host, port, certfile)

    def __setitem__(self, sym_name, input_expr):
        """register function symbol with implementation"""
        super(KamodoClient, self).__setitem__(sym_name, input_expr)
        symbol, args, lhs_units, lhs_expr = self.parse_key(sym_name)
        self._rpc_funcs[str(type(symbol))] = FunctionRPC(self[symbol], self.verbose)

    async def client_reader(self, client, reader):
        """
        Reader for the client side.
        """
        while True:
            data = await reader.read(4096)
            client.write(data)

    async def client_writer(self, client, writer):
        """
        Writer for the client side.
        """
        while True:
            data = await client.read(4096)
            writer.write(data.tobytes())
            await writer.drain()

    async def client(self, host, port, certfile):
        """
        Method to start communication as asynchronous client.
        """
        if self.verbose:
            print(f'connecting to server with {certfile}')
        try:
            ctx = ssl.create_default_context(
                ssl.Purpose.SERVER_AUTH, cafile=certfile
            )
        except FileNotFoundError:
            raise FileNotFoundError(f'{certfile} required in local directory.')

        # Handle both IPv4 and IPv6 cases
        try:
            if self.verbose:
                print("Trying IPv4")
            reader, writer = await asyncio.open_connection(
                host, port, ssl=ctx,
                family=socket.AF_INET
            )
        except Exception:
            if self.verbose:
                print("Trying IPv6")
            reader, writer = await asyncio.open_connection(
                host, port, ssl=ctx,
                family=socket.AF_INET6
            )

        if self.verbose:
            print('connection open, starting TwoPartyClient')

        # Start TwoPartyClient using TwoWayPipe (takes no arguments in this mode)
        client = capnp.TwoPartyClient()

        # Assemble reader and writer tasks, run in the background
        coroutines = [self.client_reader(client, reader), self.client_writer(client, writer)]
        asyncio.gather(*coroutines, return_exceptions=True)
        self._client = client.bootstrap().cast_as(kamodo_capnp.Kamodo)
        self._remote_fields = (await self._client.getFields().a_wait()).fields
        self._remote_math = (await self._client.getMath().a_wait()).math

        for entry in self._remote_fields.entries:
            if self.verbose:
                print('registering {}'.format(entry.key))
            await self.register_remote_field(entry)

    def connect(self, host, port, certfile):
        loop = asyncio.get_event_loop()
        loop.run_until_complete(self.client(host, port, certfile))

    async def register_remote_field(self, entry):
        """resolve the remote signature
        f(*args, **kwargs) -> f(x,y,z=value)
        """
        symbol = entry.key
        field = entry.value

        meta = field.meta
        arg_units = rpc_map_to_dict(meta.argUnits)

        defaults_ = (await field.func.getKwargs().a_wait()).kwargs
        func_defaults = {_.name: from_rpc_literal(_.value) for _ in defaults_}
        func_args_ = [str(_) for _ in (await field.func.getArgs().a_wait()).args]
        func_args = [_ for _ in func_args_ if _ not in func_defaults]

        if len(meta.equation) > 0:
            equation = meta.equation
        else:
            equation = None

        hidden_args = list(meta.hiddenArgs)

        @kamodofy(units=meta.units,
                  arg_units=arg_units,
                  citation=meta.citation,
                  equation=equation,
                  hidden_args=hidden_args)
        @forge.sign(*construct_signature(*func_args, **func_defaults))
        @wrap_async
        async def remote_func(*args, **kwargs):
            args_ = [to_rpc_literal(arg) for arg in args]
            kwargs_ = [dict(name=k, value=to_rpc_literal(v)) for k, v in kwargs.items()]
            result = (await field.func.call(args=args_, kwargs=kwargs_).a_wait()).result
            return from_rpc_literal(result)

        self[symbol] = remote_func
        self._rpc_funcs[symbol] = field.func

    async def get_remote_composition(self, expr, **kwargs):
        """Generate a callable function composition that is executed remotely"""

        async def remote_composition(**params):
            remote_expr = to_rpc_expr(expr, expressions=self._expressions, **params, **kwargs)
            evaluate_expr = self._client.evaluate(remote_expr)  # .wait()
            result_message = await evaluate_expr.value.read().a_wait()
            return from_rpc_literal(result_message.value)

        remote_composition.__name__ = str(expr)
        return remote_composition

    def vectorize_function(self, symbol, rhs_expr, composition):
        """lambdify the input expression using server-side promises"""
        if self.verbose:
            print('vectorizing {} = {}'.format(symbol, rhs_expr))
            print('composition keys {}'.format(list(composition.keys())))
        func = self.get_remote_composition(rhs_expr, **self._rpc_funcs)
        self._expressions[str(type(symbol))] = rhs_expr

        signature, defaults = sign_defaults(symbol, rhs_expr, composition)
        return signature(func)


class KamodoAPI(Kamodo):
    """JSON api wrapper for kamodo services"""

    def __init__(self, url_path, **kwargs):
        self._url_path = url_path
        super(KamodoAPI, self).__init__(**kwargs)

        self._kdata = self._get(self._url_path)

        self._defaults = {}
        self._data = {}

        self.__doc__ = requests.get('{}/doc'.format(self._url_path)).text

        for k, v in self._kdata.items():
            # get defaults for this func
            default_path = '{}/{}/defaults'.format(self._url_path, k)
            self._defaults[k] = self._get(default_path)

            # get cached data (result of calling with no args)
            data_path = '{}/{}/data'.format(self._url_path, k)
            self._data[k] = self._get(data_path)

            doc_path = '{}/{}/doc'.format(self._url_path, k)
            func_doc = requests.get(doc_path).text

            # get kamodofied function
            func = self.load_func(k)
            func.__doc__ += '\n' + func_doc
            self[k] = kamodofy(
                func,
                data=self._data[k],
                units=v['units'],
            )

    def _get(self, url_path):
        req_result = requests.get(url_path)
        try:
            result = req_result.json()  # returns a dictionary maybe
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
            params=params).json()  # returns a dictionary

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
            func = k[symbol]

            # rhs = signature['rhs']
            registry_name = '{}_{}'.format(name, kname)
            kamodo[registry_name] = func # already kamodofied

    return kamodo


def copy_func(f):
    """Based on http://stackoverflow.com/a/6528148/190597 (Glenn Maynard)
                https://stackoverflow.com/a/13503277 (unutbu)
    """
    if isinstance(f, np.vectorize):
        g = copy.deepcopy(f)
    else:
        g = types.FunctionType(f.__code__, f.__globals__, name=f.__name__,
                               argdefs=f.__defaults__,
                               closure=f.__closure__)
        g.__kwdefaults__ = f.__kwdefaults__

    g = functools.update_wrapper(g, f)
    if hasattr(f, 'meta'):
        g.meta = f.meta
    if hasattr(f, '_repr_latex_'):
        g._repr_latex_ = f._repr_latex_
    return g


def from_kamodo(kobj, **funcs):
    """copies a kamodo object, inserting additional functions"""
    knew = Kamodo()

    for name, signature in kobj.signatures.items():
        symbol = signature['symbol']
        knew[symbol] = copy_func(kobj[symbol])
    for symbol, func in funcs.items():
        knew[symbol] = func
    return knew


def get_figures(func, iterator, verbose=False):
    plots = []
    for a in func(iterator):
        if verbose:
            print('registering {}'.format(a.__name__), end=' ')
        k = Kamodo(**{a.__name__: a})
        if verbose:
            print('getting plot for {}'.format(a.__name__), end=' ')
        fig_plot = k.plot(a.__name__)
        if verbose:
            print('calling full_figure_for_development', end=' ')
        # full_fig = fig_plot.full_figure_for_development(warn=False)
        full_fig = fig_plot
        if verbose:
            print('appending {}'.format(a.__name__))
        plots.append(full_fig)
    return plots


# +
def animate(func_, iterator=None, verbose=False):
    defaults = get_defaults(func_)
    if len(defaults) > 1:
        raise NotImplementedError(
            "Animations with {} args not yet supported".format(len(defaults)))

    if iterator is None:
        param, iterator = list(defaults.items())[0]
    else:
        param = list(defaults.keys())[0]
        iterator = list(iterator)

    figures = get_figures(func_, iterator, verbose)

    print(len(figures), ' figures')

    axes_ranges = get_ranges(figures)

    layout = figures[0]['layout']
    layout.update(axes_ranges)
    # make figure

    fig_dict = {
        "data": figures[0]['data'],
        "layout": layout,
        "frames": []
    }

    fig_dict["layout"]["updatemenus"] = [
        {
            "buttons": [
                {
                    "args": [None, {"frame": {"duration": 100, "redraw": True},
                                    "fromcurrent": True,
                                    "transition": {"duration": 100,
                                                   "easing": "quadratic-in-out"}}],
                    "label": "Play",
                    "method": "animate",
                },
                {
                    "args": [[None], {"frame": {"duration": 0, "redraw": False},
                                      "mode": "immediate",
                                      "transition": {"duration": 0}}],
                    "label": "Pause",
                    "method": "animate"
                }
            ],
            "direction": "left",
            "pad": {"r": 10, "t": 87},
            "showactive": False,
            "type": "buttons",
            "x": 0.1,
            "xanchor": "right",
            "y": 0,
            "yanchor": "top"
        }
    ]

    sliders_dict = {
        "active": 0,
        "yanchor": "top",
        "xanchor": "left",
        "currentvalue": {
            "font": {"size": 20},
            "prefix": "{}: ".format(param),
            "visible": True,
            "xanchor": "right"
        },
        "transition": {"duration": 100, "easing": "cubic-in-out"},
        "pad": {"b": 10, "t": 50},
        "len": 0.9,
        "x": 0.1,
        "y": 0,
        "steps": []
    }

    # make frames
    for p, figure in zip(iterator, figures):
        frame = {"data": figure['data'],
                 "name": str(p)}
        #         print(figure['layout']['xaxis']['range'])

        fig_dict["frames"].append(frame)
        slider_step = {"args": [
            [p],
            {"frame": {"duration": 100, "redraw": True},
             "mode": "immediate",
             "transition": {"duration": 300}}],
            "label": '{:.2f}'.format(p),
            "method": "animate"}
        sliders_dict["steps"].append(slider_step)
    #     print(len(fig_dict['frames']), ' frames')

    #     fig_dict["data"] = fig_dict["data"] + fig_dict["frames"][0]["data"]
    fig_dict["layout"]["sliders"] = [sliders_dict]
    fig = go.Figure(fig_dict)
    return fig


def kamodo_constructor(loader: yaml.SafeLoader, node: yaml.nodes.MappingNode) -> Kamodo:
    """Construct a kamodo object from yaml."""
    return Kamodo(**loader.construct_mapping(node))


def kamodo_client_constructor(loader: yaml.SafeLoader, node: yaml.nodes.MappingNode) -> KamodoClient:
    """Construct an kamodo."""
    return KamodoClient(**loader.construct_mapping(node))


def yaml_loader():
    """Add Kamodo constructors to PyYAML loader."""
    loader = yaml.SafeLoader
    loader.add_constructor("!Kamodo", kamodo_constructor)
    loader.add_constructor("!KamodoClient", kamodo_client_constructor)
    return loader
