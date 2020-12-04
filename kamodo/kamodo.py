"""
Copyright © 2017 United States Government as represented by the Administrator, National Aeronautics and Space Administration.  
No Copyright is claimed in the United States under Title 17, U.S. Code.  All Other Rights Reserved.
"""
try:
    import pytest
except ImportError:
    pass

import numpy as np
from sympy import Integral, Symbol, symbols, Function

try:
    from sympy.parsing.sympy_parser import parse_expr
except ImportError:  # occurs in python3
    from sympy import sympify as parse_expr

from collections import OrderedDict
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

from sympy import Wild
from types import GeneratorType
import inspect

import re



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


# def get_expr_with_units(expr, units_map):
#     """Replaces symbols with symbol*unit"""
#     subs = []
#     for symbol, unit in list(units_map.items()):
#         if type(symbol) == str:
#             subs.append((symbol, Symbol(symbol) * unit))
#         else:  # assume symbol
#             subs.append((symbol, symbol * unit))
#     new_expr = expr.subs(subs)
#     return new_expr


# def get_expr_without_units(expr, to, units_map, dimensionless=False):
#     """Converts an expression with units to one without

#     apply conversion factors where necessary"""
#     if not dimensionless:
#         expr = sympy_units.convert_to(expr, to)
#     subs = []
#     for s in units_map:
#         if type(s) == str:
#             subs.append((Symbol(s) * to, Symbol(s)))
#         else:  # assume symbol
#             subs.append((s * to, s))
#     return expr.subs(subs)


def args_from_dict(expr, local_dict):
    """retrieve the symbols of an expression

    if a  {str: symbol} mapping is provided, use its symbols instead
    """
    args = []
    if local_dict is not None:
        for a in expr.args:
            if str(a) in local_dict:
                args.append(local_dict[str(a)])
            else:
                args.append(a)
        return tuple(args)
    else:
        return expr.args

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

    We could return symbol and dictionary mapping symbols to units

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
            f -> f, -> {f:''}

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


def parse_lhs(lhs, local_dict):
    """parse the left-hand-side

    returns: symbol, arguments, units, parsed_expr
    """
    lhs, unit_dict = extract_units(lhs)
    parsed = parse_expr(lhs)
    try:
        args = args_from_dict(parsed, local_dict)
    except:
        print(local_dict)
        raise
    symbol = expr_to_symbol(parsed, args)
    return symbol, args, unit_dict, parsed


def parse_rhs(rhs, is_latex, local_dict=None):
    """parse the right-hand-side of an equation

    subsititute variables from local_dict where available
    """
    if is_latex:
        if local_dict is not None:
            expr = parse_latex(rhs).subs(local_dict)
        else:
            expr = parse_latex(rhs)
    else:
        expr = parse_expr(rhs).subs(local_dict)
    return expr


def get_function_args(func, hidden_args=[]):
    """converts function arguments to list of symbols"""
    return symbols([a for a in getfullargspec(func).args if a not in hidden_args])


# def wildcard(expr):
#     """Replace all free symbols with the wild card symbol

#     not sure when this code is used
#     """
#     result = expr
#     for s in expr.free_symbols:
#         result = result.replace(s, Wild(str(s)))
#     return result

# def is_dimensionless(units):
#     return not hasattr(units, 'dimension')


# def validate_units(expr, units, verbose=False):
#     """Check that the right-hand-side units are consistent

#     First replace the expression units mapping.

#     IF there are free symbols remaining, assume the expression
#     is dimensionless.

#     If there are no free symbols remaining,
#     assume all symbols were replaced with their respective units
#     """

#     if hasattr(expr, 'rhs'):
#         expr_zero = expr.rhs - expr.lhs
#         result = expr_zero.subs(units, simultaneous=True)
#     else:
#         result = expr.subs(units, simultaneous=True)

#     if len(result.free_symbols) != 0:
#         return Dimension(1)

#     bases = set()
#     for symbol_ in result.free_symbols:
#         if str(symbol_) not in units:
#             raise KeyError('invalid expr {}, cannot find {} in units: {}'.format(
#                 expr, symbol_, units))
#         unit = units[str(symbol_)]
#         bases.add(unit.dimension)
#     if verbose:
#         print('bases: {}'.format(bases))
#     if len(bases) == 1:
#         print('case 1: {}'.format(expr))
#         return result
#     if len(bases) > 1:
#         print('case 2 {}'.format(expr))
#         raise ValueError('Dimension mismatch: {}, bases: {}'.format(result, bases))
#     if len(result.args) > 0:
#         print('case 3 {}'.format(expr))
#         terms = result.as_terms()[1]
#         if len(terms) > 1:
#             print('case 4 {} \n\targs: {}\n\tterms: {}'.format(expr, result.args, terms))
#             for arg_ in terms:
#                 if hasattr(arg_, 'dimension'):
#                     bases.add(arg_.dimension)
#             if len(bases) > 1:
#                 print('case 5 {}'.format(expr))
#                 raise ValueError('Dimension mismatch for {}\n\t{}\n\t{}'.format(expr, result.args, bases))
#             else:
#                 print('case 6 {}'.format(expr), bases)
#                 return expr
#         return expr

def validate_units(expr, units, verbose=False):
    """Check that the right-hand-side units are consistent

    First replace the expression units mapping.

    IF there are free symbols remaining, assume the expression
    is dimensionless.

    If there are no free symbols remaining,
    assume all symbols were replaced with their respective units
    """

    if hasattr(expr, 'rhs'):
        expr_zero = expr.rhs - expr.lhs
        result = expr_zero.subs(units, simultaneous=True)
    else:
        result = expr.subs(units, simultaneous=True)

    if len(result.free_symbols) != 0:
        return Dimension(1)

    result = result.expand()
    terms = result.as_terms()[1]
    bases = set()
    for term in terms:
        if hasattr(term, 'dimension'):
            bases.add(term)
    if len(bases) > 0:
        try:
            convert_to(result, list(bases)[0])
        except:
            raise ValueError('Cannot convert between {}'.format(bases))
    return result



# def match_dimensionless_units(lhs_units, rhs_units):
#     '''if lhs_units is dimensionless and rhs_units is not dimensionless,
#     assign rhs_units to lhs_units (and vice versa)'''
#     if lhs_units == Dimension(1):  # f = ...
#         if rhs_units != lhs_units:  # f = rho[kg/m^3]
#             lhs_units = rhs_units  # f[kg/m^3] = rho[kg/m^3]
#     elif rhs_units == Dimension(1):  # ... = x
#         if rhs_units != lhs_units:  # f[kg/m^3] = x
#             rhs_units = lhs_units  # f[kg/m^3] = x[kg/m^3]
#     return lhs_units, rhs_units


# def check_unit_compatibility(rhs_units, lhs_units):
#     """This fails to raise an error when rhs and lhs units are incompatible"""
#     try:
#         assert sympy_units.convert_to(rhs_units + lhs_units, lhs_units)
#     except:
#         raise NameError('incompatible units:{} and {}'.format(lhs_units, rhs_units))


class Kamodo(collections.OrderedDict):
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
                # elif len(components) == 3:
                #     # function(arg)[function(arg_unit) = unit] = expr
                #     if self.verbose:
                #         print('updating unit registry')
                #     lhs = self.update_unit_registry(func.strip('$'))
                #     self[lhs] = rhs

        for sym_name, expr in list(kwargs.items()):
            if self.verbose:
                print('registering {} with {}'.format(sym_name, expr))
            self[sym_name] = expr


    def register_symbol(self, symbol):
        self.symbol_registry[str(type(symbol))] = symbol

    # def load_symbol(self, sym_name):
    #     symbol = self.symbol_registry[sym_name]
    #     signature = self.signatures[str(symbol)]
    #     lhs_expr = signature['lhs']
    #     units = signature['units']
    #     return symbol, symbol.args, units, lhs_expr

    def remove_symbol(self, sym_name):
        if self.verbose:
            print('removing {} from symbol_registry'.format(sym_name))
        symbol = self.symbol_registry.pop(sym_name)
        if self.verbose:
            print('removing {} from signatures'.format(symbol))
        self.signatures.pop(str(symbol))
        if self.verbose:
            print('removing {} from unit_registry'.format(symbol))
        remove = []
        for sym in self.unit_registry:
            if is_function(sym): # {rho(x): rho(cm), rho(cm): kg}
                if type(sym) == type(symbol):
                    remove.append(sym)
        for sym in remove:
            self.unit_registry.pop(sym)

        # if sym_name in self.unit_registry:
        #     func_unit = self.unit_registry.pop(sym_name) # rho(x): rho(cm)
        #     if func_unit in self.unit_registry:
        #         self.unit_registry.pop(func_unit) # rho(cm): kg/m^3

    def parse_key(self, sym_name):
        args = tuple()
        if sym_name not in self:
            if self.verbose:
                print('parsing new key {}'.format(sym_name))
            if type(sym_name) is str:
                sym_name = sym_name.strip('$').strip()
                if sym_name not in self.symbol_registry:
                    symbol, args, unit_dict, lhs_expr = parse_lhs(sym_name, self.symbol_registry)
                    symbol_str = str(symbol).replace(' ', '')
                    if symbol_str in unit_dict:
                        units = unit_dict[symbol_str]
                    elif str(type(symbol)) in unit_dict:
                        units = unit_dict[str(type(symbol))]
                    else:
                        raise NameError('{} not found in unit_dict {}'.format(symbol_str, unit_dict))
                    if self.verbose:
                        print('newly parsed symbol:', symbol, type(symbol))
                    if str(type(symbol)) in self.symbol_registry:
                        raise KeyError("{} found in symbol_registry".format(str(type(symbol))))
                else:
                    raise KeyError("{} found in symbol_registry".format(sym_name))
            else:
                if type(sym_name) is Symbol:
                    symbol = sym_name
                    units = ''  # where else do we get these?
                    lhs_expr = symbol
                else:
                    # cast the lhs into a string and parse it
                    symbol, args, unit_dict, lhs_expr = parse_lhs(str(sym_name), self.symbol_registry)
                    symbol_str = str(symbol).replace(' ', '')
                    if symbol_str in unit_dict:
                        units = unit_dict[symbol_str]
                    elif str(type(symbol)) in unit_dict:
                        units = unit_dict[str(type(symbol))]
                    else:
                        raise NameError('{} not found in unit_dict {}'.format(symbol_str, unit_dict))
        else:
            symbol = sym_name
            if self.verbose:
                print('{} found in keys'.format(sym_name))
            try:
                units = self.signatures[sym_name]['units']
                lhs_expr = self.signaturs[sym_name]['lhs']
            except KeyError:
                units = ''
                lhs_expr = symbol
        return symbol, args, units, lhs_expr

    def parse_value(self, rhs_expr, local_dict=None):
        """returns an expression from string"""
        if type(rhs_expr) is str:
            is_latex = '$' in rhs_expr
            rhs_expr = rhs_expr.strip('$').strip()
            rhs_expr = parse_rhs(rhs_expr, is_latex, local_dict=local_dict)
        return rhs_expr

    def check_or_replace_symbol(self, symbol, free_symbols, rhs_expr=None):
        """Rules of replacement:

        """
        lhs_args = set(symbol.args)
        if lhs_args != set(free_symbols):
            free_symbols_ = tuple(free_symbols)
            if len(free_symbols) == 1:
                try:
                    symbol = parse_expr(str(type(symbol)) + str(free_symbols_))
                except:
                    raise NotImplementedError('cannot parse', str(type(symbol)) + str(free_symbols_))
            else:
                if len(lhs_args) > 0:
                    raise NameError("Mismatched symbols {} and {}".format(lhs_args, free_symbols))
                else:
                    if self.verbose:
                        print('type of rhs symbols:', type(free_symbols))
                    if type(free_symbols) == set:
                        free_symbols_ = sort_symbols(free_symbols)
                    if self.verbose:
                        print('replacing {} with {}'.format(
                            symbol, str(type(symbol)) + str(free_symbols_)))
                    try:
                        symbol = parse_expr(str(type(symbol)) + str(free_symbols_))
                    except:
                        raise NotImplementedError('cannot parse', str(type(symbol)) + str(free_symbols_))

        return symbol

    def validate_function(self, lhs_expr, rhs_expr):
        assert lhs_expr.free_symbols == rhs_expr.free_symbols

    def get_composition(self, lhs_expr, rhs_expr):
        composition = dict()
        for k in list(self.keys()):
            if len(rhs_expr.find(k)) > 0:
                if self.verbose:
                    print('composition detected: found {} in {} = {}'.format(k, lhs_expr, rhs_expr))
                composition[str(k)] = self[k]
            else:
                if self.verbose:
                    print('{} {} not in {} = {}'.format(k, type(k), lhs_expr, rhs_expr))
        return composition

    def vectorize_function(self, symbol, rhs_expr, composition):
        try:
            func = lambdify(symbol.args, rhs_expr, modules=['numexpr'])
            if self.verbose:
                print('lambda {} = {} labmdified with numexpr'.format(symbol.args, rhs_expr))
        except:
            func = lambdify(symbol.args, rhs_expr, modules=['numpy', composition])
            if self.verbose:
                print('lambda {} = {} lambdified with numpy and composition:'.format(
                    symbol.args, rhs_expr))
                for k, v in composition.items():
                    print('\t', k, v)
        return func

    def update_unit_registry(self, func, arg_units=None):
        """Inserts unit functions into registry"""
        lhs, unit_dict = extract_units(func)
        if arg_units is None:
            arg_units = {}
        for key, value in unit_dict.items():
            if key != lhs:
                arg_units[parse_expr(key)] = get_unit(value)

        lhs_expr = parse_expr(lhs)
        func_units = lhs_expr.subs(arg_units)
        self.unit_registry[lhs_expr] = func_units
        output_units = unit_dict[lhs]
        self.unit_registry[func_units] = get_unit(output_units)
        return lhs

        # func_symbol = func.split('[')[0]
        # unit = func.split('[')[1].split(']')[0]
        # rhs = func.split('=')[-1]
        # if '=' in unit:
        #     unit_lhs, unit_rhs = unit.split('=')
        #     unit_lhs = parse_expr(unit_lhs).subs(unit_subs)
        #     self.unit_registry[parse_expr(func_symbol)] = unit_lhs
        #     unit_rhs = get_unit(unit_rhs)
        #     self.unit_registry[unit_lhs] = unit_rhs
        # else:
        #     self.unit_registry[parse_expr(func_symbol)] = get_unit(unit)
        # return func_symbol, rhs

    def register_signature(self, symbol, units, lhs_expr, rhs_expr):
        if isinstance(units, str):
            unit_str = units
            if self.verbose:
                print('unit str {}'.format(unit_str))
        else:
            if self.verbose:
                print('getting abbreviation for', units)
            unit_abbrev = get_abbrev(units)
            if unit_abbrev is not None:
                unit_str = str(unit_abbrev)
            else:
                unit_str = None
        self.signatures[str(symbol)] = dict(
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
            except:  # will not work on methods
                pass

        arg_units = func.meta.get('arg_units', None)
        if arg_units is not None:
            unit_dict = {arg: get_unit(arg_unit) for arg, arg_unit in arg_units.items()}
            unit_expr = lhs_symbol.subs(unit_dict)
            self.unit_registry[lhs_symbol] = unit_expr
            self.unit_registry[unit_expr] = get_unit(units)
        else:
            self.unit_registry[lhs_symbol] = get_unit(units)
        super(Kamodo, self).__setitem__(lhs_symbol, func)  # assign key 'f(x)'
        self.register_signature(lhs_symbol, units, lhs_expr, rhs)
        super(Kamodo, self).__setitem__(type(lhs_symbol), self[lhs_symbol])  # assign key 'f'
        self.register_symbol(lhs_symbol)

    # def check_consistency(self, input_expr, units):
    #     # check that rhs units are consistent
    #     rhs_expr = self.parse_value(input_expr, self.symbol_registry)
    #     units_map = self.get_units_map()
    #     lhs_units = get_unit(units)
    #     rhs_units = validate_units(rhs_expr, units_map, self.verbose)
    #     if self.verbose:
    #         print('rhs_units: {}'.format(rhs_units))
    #     try:
    #         lhs_units, rhs_units = match_dimensionless_units(lhs_units, rhs_units)
    #         check_unit_compatibility(rhs_units, lhs_units)
    #     except:
    #         print(type(rhs_units))
    #         print(get_unit(units), validate_units(rhs_expr, units_map, self.verbose))
    #         raise

    #     rhs_expr_with_units = get_expr_with_units(rhs_expr, units_map)

    #     if self.verbose:
    #         print('rhs_expr with units:', rhs_expr_with_units)

    #     # convert back to expression without units for lambdify
    #     if units != '':
    #         try:
    #             if self.verbose:
    #                 print('converting to {}'.format(lhs_units))
    #                 for k, v in list(units_map.items()):
    #                     print('\t', k, v, type(k))
    #             rhs_expr = get_expr_without_units(
    #                 rhs_expr_with_units,
    #                 lhs_units,
    #                 units_map,
    #                 dimensionless=is_dimensionless(rhs_units))
    #             if self.verbose:
    #                 print('rhs_expr without units:', rhs_expr)
    #         except:
    #             print('error with units? [{}]'.format(units))
    #             raise
    #     else:
    #         if lhs_units != Dimension(1):  # lhs_units were obtained from rhs_units
    #             units = str(lhs_units)
    #     return units, rhs_expr

    def __setitem__(self, sym_name, input_expr):
        """Assigns a function or expression to a new symbol,
        performs unit conversion where appropriate

        """
        if not isinstance(sym_name, str):
            sym_name = str(sym_name)

        if self.verbose:
            print('')
        try:
            symbol, args, lhs_units, lhs_expr = self.parse_key(sym_name)
        except KeyError as error:
            if self.verbose:
                print('could not use parse_key with {}'.format(sym_name))
                print(error)
            found_sym_name = str(error).split('found')[0].strip("'").strip(' ')
            if self.verbose:
                print('replacing {}'.format(found_sym_name))
            self.remove_symbol(found_sym_name)
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
                if symbol in self.unit_registry:
                    units = get_expr_unit(symbol, self.unit_registry)
                    if self.verbose:
                        print('{} has units {}'.format(sym_name, units))
                else:
                    if self.verbose:
                        print(sym_name,
                              symbol,
                              'had no units. Getting units from {}'.format(rhs_expr))

                    expr_unit = get_expr_unit(rhs_expr, self.unit_registry, self.verbose)
                    arg_units = get_arg_units(rhs_expr, self.unit_registry)

                    if expr_unit == Dimension(1):
                        expr_unit = None

                    if self.verbose:
                        print('registering {} with {} {}'.format(symbol, expr_unit, arg_units))

                    if (symbol not in self.unit_registry) and (expr_unit is not None):
                        self.unit_registry[symbol] = symbol.subs(arg_units)
                        self.unit_registry[symbol.subs(arg_units)] = expr_unit

                    if is_function(expr_unit):
                        self.unit_registry[expr_unit] = get_expr_unit(
                            expr_unit,
                            self.unit_registry,
                            self.verbose)

                    if expr_unit is not None:
                        lhs_units = str(get_abbrev(get_expr_unit(expr_unit, self.unit_registry, self.verbose)))

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
                    print('\t{}: {}'.format(k, v))

            units = get_expr_unit(symbol, self.unit_registry)
            units = get_abbrev(units)
            if units is not None:
                units = str(units)
            else:
                units = ''
            if self.verbose:
                print('units after resolve', symbol, units)
                for k, v in self.unit_registry.items():
                    print('\t{}: {}'.format(k, v))

            rhs_args = rhs_expr.free_symbols

            try:
                symbol = self.check_or_replace_symbol(symbol, rhs_args, rhs_expr)
                self.validate_function(symbol, rhs_expr)

            except:
                if self.verbose:
                    print('\n Error in __setitem__', input_expr)
                    print(symbol, lhs_expr, rhs_args)
                    print('symbol registry:', self.symbol_registry)
                    print('signatures:', self.signatures)
                    print('unit registry:', self.unit_registry)
                raise

            # composition = self.get_composition(lhs_expr, rhs_expr)
            composition = {str(k_): self[k_] for k_ in self}
            arg_units = {}
            if symbol in self.unit_registry:
                unit_args = self.unit_registry[symbol]
                if unit_args is not None:
                    if len(unit_args.args) == len(symbol.args):
                        for arg, unit in zip(symbol.args, unit_args.args):
                            arg_units[str(arg)] = str(get_abbrev(unit))
            func = kamodofy(
                self.vectorize_function(symbol, rhs_expr, composition),
                units=units,
                arg_units=arg_units)
            self.register_signature(symbol, units, lhs_expr, rhs_expr)
            super(Kamodo, self).__setitem__(symbol, func)
            super(Kamodo, self).__setitem__(type(symbol), self[symbol])
            self.register_symbol(symbol)
            # self[symbol].meta = dict(units=units)


    def __getitem__(self, key):
        try:
            return super(Kamodo, self).__getitem__(key)
        except:
            return self[self.symbol_registry[key]]

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            name_ = symbols(name, cls=Function)
            if name_ in self:
                return self[name_]
        raise AttributeError(name)

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such field: " + name)

    def get_units_map(self):
        """Maps from string units to symbolic units"""
        d = dict()
        star = Wild('star')
        for k, v in list(self.signatures.items()):
            unit = get_unit(v['units'])
            d[k] = unit
            d[v['symbol']] = unit

        return d

    def to_latex(self, keys=None, mode='equation'):
        """Generate list of LaTeX-formated formulas"""
        if keys is None:
            keys = list(self.signatures.keys())
        repr_latex = ""
        for k in keys:
            lhs = self.signatures[k]['symbol']
            rhs = self.signatures[k]['rhs']
            units = self.signatures[k]['units']
            arg_units = self.unit_registry.get(lhs, None)

            # f(x[cm],y[km],z)[km] = expr
            units = self.unit_registry.get(arg_units)
            if units is not None:
                units = '{}'.format(get_abbrev(units))
            else:
                units = ''

            if arg_units is not None:
                lhs_str = "{}({})".format(
                    latex(type(lhs)),
                    ",".join([
                        "{}[{}]".format(latex(arg), get_abbrev(unit)) for arg, unit in zip(lhs.args, arg_units.args)
                        ])
                    )
            else:
                lhs_str = latex(lhs)

            if len(units) > 0:
                lhs_str += "[{}]".format(units)


            if type(rhs) == str:
                latex_eq = rhs
                # latex_eq = latex(Eq(lhs, parse_latex(rhs)), mode = mode)
            else:
                if rhs is not None:
                    try:
                        latex_eq = latex(Eq(lhs, rhs), mode=mode)
                        latex_eq_rhs = latex(rhs) # no $$ delimiter
                    except:
                        lambda_ = symbols('lambda', cls=UndefinedFunction)
                        latex_eq = latex(Eq(lhs, lambda_(*lhs.args)), mode=mode)
                        latex_eq_rhs = latex(lambda_(*lhs.args)) # no $$
                else:
                    lambda_ = symbols('lambda', cls=UndefinedFunction)
                    latex_eq = latex(Eq(lhs, lambda_(*lhs.args)), mode=mode)
                    latex_eq_rhs = latex(lambda_(*lhs.args)) # no $$
            if units is None:
                units = ''
            if len(str(units)) > 0:
                # latex_eq = latex_eq.replace('=','\\text{' + '[{}]'.format(units) + '} =')
                latex_eq = latex_eq.replace('=', '[{}] ='.format(units))

            repr_latex += r"\begin{equation}"
            repr_latex += "{} = {}".format(lhs_str, latex_eq_rhs)
            repr_latex += r"\end{equation}"

            # repr_latex += latex_eq

        return beautify_latex(repr_latex).encode('utf-8').decode()

    def _repr_latex_(self):
        """Provide notebook rendering of formulas"""
        return self.to_latex()

    def detail(self):
        """Constructs a pandas dataframe from signatures"""
        return pd.DataFrame(self.signatures).T

    def get_signature(self, name):
        """Get the signature for the named variable"""
        return self.signatures[str(self.symbol_registry[name])]

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
        if params is not None:
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
        result = self.evaluate(variable, **kwargs)
        signature = self.get_signature(variable)
        units = signature['units']
        if units != '':
            units = '[{}]'.format(units)
        title = self.to_latex([str(self.symbol_registry[variable])], mode='inline')  # .replace('\\operatorname','')
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
        try:
            out_dim, arg_dims = get_plot_key(result[variable].shape, *arg_shapes)
        except:
            print(arg_shapes)
            raise
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


def compose(**kamodos):
    """Kamposes multiple kamodo instances into one"""
    kamodo = Kamodo()
    for kname, k in kamodos.items():
        for name, symbol in k.symbol_registry.items():
            signature = k.signatures[str(symbol)]
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
    for var_symbol, func in kobj.items():
        if not is_function(var_symbol):
            knew[var_symbol] = func
    for var_symbol, func in funcs.items():
        knew[var_symbol] = func
    return knew
