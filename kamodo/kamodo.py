"""
Copyright © 2017 United States Government as represented by the Administrator, National Aeronautics and Space Administration.  
No Copyright is claimed in the United States under Title 17, U.S. Code.  All Other Rights Reserved.
"""
try:
	import pytest
except ImportError:
	pass

from collections import OrderedDict
import collections

import numpy as np
from sympy import Integral, Symbol, symbols, Function

from sympy.parsing.sympy_parser import parse_expr

from sympy import lambdify
from sympy.parsing.latex import parse_latex
from sympy import latex
from sympy.core.function import UndefinedFunction
from inspect import getargspec
from sympy import Eq
import pandas as pd
from IPython.display import Latex


from sympy.physics import units as sympy_units
from sympy.physics.units import Quantity
from sympy.physics.units import Dimension
from sympy import Expr
from sympy import sympify

import functools
from .util import kamodofy
from .util import sort_symbols
from .util import simulate
from .util import unit_subs
from .util import get_defaults, valid_args, eval_func
from .util import to_arrays, cast_0_dim
from .util import beautify_latex, arg_to_latex
from .util import concat_solution

import sympy.physics.units as u


import plotly.graph_objs as go
from plotly import figure_factory as ff

from plotting import *
from .util import existing_plot_types

from sympy import Wild
from types import GeneratorType

def get_unit_quantities():
	subs = {} 
	for k, v in list(sympy_units.__dict__.items()): 
		if isinstance(v, Expr) and v.has(sympy_units.Unit): 
			subs[Symbol(k)] = v # Map the `Symbol` for a unit to the quantity
	return subs


def clean_unit(unit_str):
	"""remove brackets and all white spaces in and around unit string"""
	return unit_str.strip().strip('[]').strip()

def get_unit(unit_str, unit_subs = unit_subs):
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
		raise ValueError("Unsupported unit: {} {}".format(unit_str, type(unit_str)))
	return unit

def get_expr_with_units(expr, units_map):
	"""Replaces symbols with symbol*unit"""
	subs = []
	for s,u in list(units_map.items()):
		if type(s) == str:
			subs.append((s, Symbol(s)*u))
		else: #assume symbol
			subs.append((s,s*u))
	new_expr = expr.subs(subs)
	return new_expr

def get_expr_without_units(expr, to, units_map):
	"""Converts an expression with units to one without, applying conversion factors where necessary"""
	expr = sympy_units.convert_to(expr,to)
	subs = []
	for s in units_map:
		if type(s) == str:
			subs.append((Symbol(s)*to, Symbol(s)))
		else: #assume symbol
			subs.append((s*to, s))
	return expr.subs(subs)


def args_from_dict(expr, local_dict):
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

def extract_units(lhs):
	if '[' in lhs:
		lhs, units = lhs.split('[')
		return lhs, units.split(']')[0]
	else:
		return lhs, ''

def expr_to_symbol(expr, args):
	"""f,args -> f(args) of class function, f(x) -> f(x)"""
	if type(expr) == Symbol:
		return parse_expr(str(expr) + str(args))
		# return symbols(str(expr), cls = Function)
	else:
		return expr 

def parse_lhs(lhs, local_dict):	
	lhs, units = extract_units(lhs)
	parsed = parse_expr(lhs)
	try:
		args = args_from_dict(parsed, local_dict)
	except:
		print(local_dict)
		raise
	symbol = expr_to_symbol(parsed, args)
	return symbol, args, units, parsed

def parse_rhs(rhs, is_latex, local_dict = None):
	parsed = dict()
	if is_latex:
		if local_dict is not None:
			expr = parse_latex(rhs).subs(local_dict)
		else:
			expr = parse_latex(rhs)
	else:
		expr = parse_expr(rhs).subs(local_dict)
	return expr

def parse_function(function, local_dict = None):
	parsed = dict()

	is_latex = False

	if "$" in function:
		is_latex = True
		function = function.strip("$").strip()
	
	split = [a.strip() for a in function.split('=')]
	
	if len(split) == 2:
		lhs, rhs = split
		expr = parse_rhs(rhs, is_latex, local_dict)
	else:
		lhs = split[0]
		expr = None

	symbol, args, units, lhs_expr = parse_lhs(lhs, local_dict)
	if expr is not None:
		if set(args) != expr.free_symbols:
			raise NameError("Mismatched symbols: {} and {}".format(args, expr.free_symbols))

	parsed['expr'] = expr
	parsed['symbol'] = symbol
	parsed['args'] = args
	parsed['units'] = units

	return parsed

def get_function_args(func, hidden_args = []):
	"""converts function arguments to list of symbols"""
	return symbols([a for a in getargspec(func).args if a not in hidden_args])

def wildcard(expr):
	result = expr
	for s in expr.free_symbols:
		result = result.replace(s, Wild(str(s)))
	return result


def validate_units(expr, units):
	result = expr.subs(units, simultaneous = True)
	if len(result.free_symbols) != 0:
		return Dimension(1)
	else:
		for s, unit in list(units.items()):
			try:
				result = result.replace(wildcard(s), unit)
			except:
				pass
		return result




def match_dimensionless_units(lhs_units, rhs_units):
	'''if lhs_units is dimensionless and rhs_units is not dimensionless, assign rhs_units to lhs_units (and vice versa)'''
	if lhs_units == Dimension(1): # f = ...
		if rhs_units != lhs_units: # f = rho[kg/m^3]
			lhs_units = rhs_units # f[kg/m^3] = rho[kg/m^3]
	elif rhs_units == Dimension(1): # ... = x
		if rhs_units != lhs_units: # f[kg/m^3] = x
			rhs_units = lhs_units # f[kg/m^3] = x[kg/m^3]
	return lhs_units, rhs_units
	

def check_unit_compatibility(rhs_units, lhs_units):
	try:
		assert sympy_units.convert_to(rhs_units + lhs_units, lhs_units)
	except:
		raise NameError('incompatible units:{} and {}'.format(lhs_units, rhs_units))


class Kamodo(collections.OrderedDict):
	'''Kamodo base class demonstrating common API for space weather models


	This API provides access to space weather fields and their properties through:

		interpolation of variables at user-defined points
		unit conversions
		coordinate transformations specific to space weather domains

	Required methods that have not been implemented in child classes will raise a NotImplementedError
	'''
	def __init__(self, *funcs, **kwargs):
		"""Base initialization method

		Args:
			param1 (str, optional): Filename of datafile to interpolate from

		"""

		super(Kamodo, self).__init__()
		self.symbol_registry = OrderedDict()

		symbol_dict = kwargs.pop('symbol_dict', None)

		self.verbose = kwargs.pop('verbose', False)
		self.signatures = OrderedDict()

		for func in funcs:
			if type(func) is str:
				lhs,rhs = func.split('=')
				self[lhs.strip('$')] = rhs

		for sym_name, expr in list(kwargs.items()):
			self[sym_name] = expr


	def __getitem__(self, key):
		try:
			return super(Kamodo, self).__getitem__(key)
		except:
			return self[self.symbol_registry[key]]

	def register_symbol(self, symbol):
		self.symbol_registry[str(type(symbol))] = symbol

	def load_symbol(self, sym_name):
		symbol = self.symbol_registry[sym_name]
		signature = self.signatures[str(symbol)]
		lhs_expr = signature['lhs']
		units = signature['units']
		return symbol, symbol.args, units, lhs_expr

	def remove_symbol(self, sym_name):
		if self.verbose:
			print('removing {} from symbol_registry'.format(sym_name))
		symbol = self.symbol_registry.pop(sym_name)
		if self.verbose:
			print('removing {} from signatures'.format(symbol))
		self.signatures.pop(str(symbol))

	def parse_key(self, sym_name):
		args = tuple()
		if sym_name not in self:
			if self.verbose:
				print('parsing new key {}'.format(sym_name))
			if type(sym_name) is str:
				sym_name = sym_name.strip('$').strip()
				if sym_name not in self.symbol_registry:
					symbol, args, units, lhs_expr = parse_lhs(sym_name, self.symbol_registry)
					if self.verbose:
						print('newly parsed symbol:', symbol, type(symbol))
					if str(type(symbol)) in self.symbol_registry:
						raise KeyError("{} found in symbol_registry".format(str(type(symbol))))
				else:
					raise KeyError("{} found in symbol_registry".format(sym_name))
			else:
				if type(sym_name) is Symbol:
					symbol = sym_name
					units = '' # where else do we get these?
					lhs_expr = symbol
				else:
					#cast the lhs into a string and parse it
					symbol, args, units, lhs_expr = parse_lhs(str(sym_name), self.symbol_registry) 
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


	def parse_value(self, rhs_expr, local_dict = None):
		if type(rhs_expr) is str:
			is_latex = '$' in rhs_expr
			rhs_expr = rhs_expr.strip('$').strip()
			rhs_expr = parse_rhs(rhs_expr, is_latex, local_dict = local_dict)
		return rhs_expr


	def check_or_replace_symbol(self, symbol, free_symbols, rhs_expr = None):
		"""Rules of replacement:

		"""
		lhs_args = set(symbol.args)
		if lhs_args != set(free_symbols):
			free_symbols_ = tuple(free_symbols)
			if len(free_symbols) == 1:
				symbol = parse_expr(str(type(symbol)) + str(free_symbols_))
			else:
				if len(lhs_args) > 0:
					raise NameError("Mismatched symbols {} and {}".format(lhs_args, free_symbols))
				else:
					if self.verbose:
						print('type of rhs symbols:', type(free_symbols))
					if type(free_symbols) == set:
						free_symbols_ = sort_symbols(free_symbols)
					if self.verbose:
						print('replacing {} with {}'.format(symbol, str(type(symbol)) + str(free_symbols_)))
					symbol = parse_expr(str(type(symbol)) + str(free_symbols_))

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
		return composition


	def vectorize_function(self, symbol, rhs_expr, composition):
		try:
			func = lambdify(symbol.args, rhs_expr, modules = ['numexpr'])
			if self.verbose:
				print('lambda {} = {} labmdified with numexpr'.format(symbol.args, rhs_expr))
		except:
			func = lambdify(symbol.args, rhs_expr, modules = ['numpy', composition])
			if self.verbose:
				print('lambda {} = {} lambdified with numpy and {}'.format(symbol.args, rhs_expr, composition))
		return func

	def register_signature(self, symbol, units, lhs_expr, rhs_expr):
		self.signatures[str(symbol)] = dict(
				symbol = symbol, 
				units = units, 
				lhs = lhs_expr,
				rhs = rhs_expr,
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
		if str(rhs_args[0]) == 'self': # in case function is a class method
			rhs_args.pop(0)
		
		lhs_symbol = self.check_or_replace_symbol(lhs_symbol, rhs_args)
		units = lhs_units
		if hasattr(func, 'meta'):
			if len(lhs_units) > 0:
				if lhs_units != func.meta['units']:
					raise NameError('Units mismatch:{} != {}'.format(lhs_units, func.meta['units']))
			else:
				units = func.meta['units']
			rhs = func.meta.get('equation', func)
		else:
			rhs = func
			try:
				setattr(func, 'meta', dict(units = lhs_units))
			except: # will not work on methods
				pass
		super(Kamodo, self).__setitem__(lhs_symbol, func) # assign key 'f(x)'
		self.register_signature(lhs_symbol, units, lhs_expr, rhs)
		super(Kamodo, self).__setitem__(type(lhs_symbol), self[lhs_symbol]) # assign key 'f'
		self.register_symbol(lhs_symbol) 


	def __setitem__(self, sym_name, input_expr):
		"""Assigns a function or expression to a new symbol, performing unit conversion where necessary

		"""
		if self.verbose:
			print('')
		try:
			symbol, args, units, lhs_expr = self.parse_key(sym_name)
		except KeyError as error:
			if self.verbose:
				print(error)
			found_sym_name = str(error).split('found')[0].strip("'").strip(' ')
			if self.verbose:
				print('replacing {}'.format(found_sym_name))
			self.remove_symbol(found_sym_name)
			symbol, args, units, lhs_expr = self.parse_key(sym_name)

		if hasattr(input_expr, '__call__'):
			self.register_function(input_expr, symbol, lhs_expr, units)
		else:
			rhs_expr = self.parse_value(input_expr, self.symbol_registry)
			units_map = self.get_units_map()
			lhs_units = get_unit(units)
			rhs_units = validate_units(rhs_expr, units_map) # check that rhs units are consistent

			try:
				lhs_units, rhs_units = match_dimensionless_units(lhs_units, rhs_units)
				check_unit_compatibility(rhs_units, lhs_units)
			except:
				print(type(rhs_units))
				print(get_unit(units), validate_units(rhs_expr, units_map))
				raise

			rhs_expr_with_units = get_expr_with_units(rhs_expr, units_map)

			if self.verbose:
				print('rhs_expr with units:', rhs_expr_with_units)

			# convert back to expression without units for lambdify
			if units != '':
				try:
					if self.verbose:
						print('converting to {}'.format(lhs_units))
						for k,v in list(units_map.items()):
							print('\t', k,v, type(k))
					rhs_expr = get_expr_without_units(rhs_expr_with_units, lhs_units, units_map)
					if self.verbose:
						print('rhs_expr without units:', rhs_expr)
				except:
					print('error with units? [{}]'.format(units))
					raise
			else:
				if lhs_units != Dimension(1): # lhs_units were obtained from rhs_units
					units = str(lhs_units)

			rhs_args = rhs_expr.free_symbols

			if self.verbose:
				print('lhs_expr:', lhs_expr, 'units:', lhs_units)
				print('rhs_expr:', rhs_expr, 'units:', rhs_units)
			try:
				symbol = self.check_or_replace_symbol(symbol, rhs_args, rhs_expr)
				self.validate_function(symbol, rhs_expr)

			except:
				if self.verbose:
					print('\n Error in __setitem__', input_expr)
					print(symbol, lhs_expr, rhs_args)
					print('symbol registry:', self.symbol_registry)
					print('signatures:', self.signatures)
				raise
		

			composition = self.get_composition(lhs_expr, rhs_expr)
			func = self.vectorize_function(symbol, rhs_expr, composition)
			self.register_signature(symbol, units, lhs_expr, rhs_expr)
			super(Kamodo, self).__setitem__(symbol, func)
			super(Kamodo, self).__setitem__(type(symbol), self[symbol])
			self.register_symbol(symbol)
			self[symbol].meta = dict(units = units)

	def __getattr__(self, name):
		name_ = symbols(name, cls = Function)
		if name_ in self:
			return self[name_]
		else:
			raise AttributeError("No such symbol: {}".format(name))

	def __delattr__(self, name):
		if name in self:
			del self[name]
		else:
			raise AttributeError("No such field: " + name)

	def get_units_map(self):
		"""Maps from string units to symbolic units"""
		d = dict()
		star = Wild('star')
		for k,v in list(self.signatures.items()):
			unit = get_unit(v['units'])
			d[k] = unit
			d[v['symbol']] = unit

		return d 


	def to_latex(self, keys = None, mode = 'equation'):
		"""Generate list of LaTeX-formated formulas"""
		if keys is None:
			keys = list(self.signatures.keys())
		repr_latex = ""
		for k in keys:
			lhs = self.signatures[k]['symbol']
			rhs = self.signatures[k]['rhs']
			units = self.signatures[k]['units']
			if type(rhs) == str:
				latex_eq = rhs
				# latex_eq = latex(Eq(lhs, parse_latex(rhs)), mode = mode)
			else:
				try:
					latex_eq = latex(Eq(lhs, rhs), mode = mode)
				except:
					lambda_ = symbols('lambda', cls = UndefinedFunction)
					latex_eq = latex(Eq(lhs, lambda_(*lhs.args)), mode = mode)
			if len(units) > 0:
				# latex_eq = latex_eq.replace('=','\\text{' + '[{}]'.format(units) + '} =')
				latex_eq = latex_eq.replace('=','[{}] ='.format(units))
			repr_latex += latex_eq

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
		if type(self[variable]) is np.vectorize:
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
			dense_output = True, # generate a callable solution
			events = None, # stop when event is triggered
			vectorized = True,):
		from scipy.integrate import solve_ivp
		result = solve_ivp(self[fprime], interval, y0, 
			dense_output = dense_output, 
			events = events, 
			vectorized = vectorized)
		if self.verbose:
			print(result['message'])
		def solution(t = result['t']):
			return result['sol'].__call__(t).T.squeeze()
		return solution


	def figure(self, variable, indexing = 'ij', return_type = False, **kwargs):
		result = self.evaluate(variable, **kwargs)
		signature = self.get_signature(variable)
		units = signature['units']
		if units != '':
			units = '[{}]'.format(units)
		title = self.to_latex([str(self.symbol_registry[variable])], mode = 'inline') #.replace('\\operatorname','')
		title_lhs = title.split(' =')[0]+'$'
		title_short = '{}'.format(variable + units) # something wrong with colorbar latex being vertical
		titles = dict(title = title, title_lhs = title_lhs, title_short = title_short, units = units, variable = variable)
		fig = dict()
		chart_type = None
		traces = []

		arg_shapes = get_arg_shapes(*[v for v in list(result.values())[:-1]])
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
		traces, chart_type, layout = plot_func(result, titles,indexing = indexing, verbose = self.verbose)
		# if len(result) == 2: 
		# 	arg0, val0 = result.items()[0]
		# 	if self.verbose:
		# 		print '1 input: {}-d shape: {}'.format(len(val0.shape), val0.shape)
		# 	if len(val0.shape) == 1:
		# 		lines, chart_type, layout = line_plot(result, titles, self.verbose)
		# 		traces.extend(lines)

		# 	elif len(val0.shape) == 2:
		# 		vectors, chart_type, layout = vector_plot(result, titles, self.verbose)
		# 		traces.extend(vectors)
		# 	else: 
		# 		raise NotImplementedError('shape not supported'.format(val0.shape))
		# elif len(result) == 3: 
		# 	if len(result[variable].shape) == 2:
		# 		contours, chart_type, layout = contour_plot(result, titles, indexing, self.verbose)
		# 		traces.extend(contours)

			
		# elif len(result) == 4:
		# 	result = to_arrays(result) # loses the date-time index
		# 	arg0, val0 = result.items()[0]
		# 	arg1, val1 = result.items()[1]
		# 	arg2, val2 = result.items()[2]
		# 	if self.verbose:
		# 		print '3 inputs: {}{}, {}{}, {}{}'.format(
		# 			arg0, val0.shape, 
		# 			arg1, val1.shape,
		# 			arg2, val2.shape)
		# 	if len(result[variable].shape) == 3:
		# 		planes, chart_type, layout = plane(result, titles, indexing, self.verbose)
		# 		traces.extend(planes)

		# 	elif len(result[variable].shape) == 2:
		# 		if self.verbose:
		# 			print '\t2-d output', result[variable].shape
		# 		surfaces, chart_type, layout = surface(result, titles, self.verbose)
		# 		traces.extend(surfaces)

		# 	elif len(result[variable].shape) == 1:
		# 		if self.verbose:
		# 			print '\tscalar output', result[variable].shape
		# 		if result[variable].shape == val0.shape:
		# 			lines, chart_type, layout = line_plot(result, titles, self.verbose)
		# 			traces.extend(lines)

		# 		elif result[variable].shape[0] == 1:
		# 			surfaces, chart_type, layout = surface(result, titles, self.verbose)
		# 			traces.extend(surfaces)

		# 	elif len(result[variable]) == 1: #scalar or color value?
		# 		surfaces, chart_type, layout = surface(result, titles, self.verbose)
		# 		traces.extend(surfaces)

		# 	else:
		# 		print '{}-d not yet supported'.format(result[variable].shape)
			
		# else:
		# 	print 'not supported yet: {}'.format(type(result[variable]))


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
			fig = self.figure(variable, return_type = True, **kwargs)
			if fig['chart_type'] is None:
				raise AttributeError("No chart_type for this trace")
			else:
				if self.verbose:
					print('chart type:', fig['chart_type'])
				return go.Figure(data = fig['data'], layout = fig['layout'])
		else:
			traces = []
			layouts = []
			for variable, kwargs in list(figures.items()):
				fig = self.figure(variable, **kwargs)
				traces.extend(fig['data'])
				layouts.append(fig['layout'])
			return go.Figure(data = traces, layout = layouts[-1])


##------------ Tests Below --------------------

def test_Kamodo_expr():
	a, b, c, x, y, z = symbols('a b c x y z')
	kamodo = Kamodo(a = x**2, verbose = False)
	try:
		assert kamodo.a(3) == 3**2
	except:
		print(kamodo.symbol_registry)
		raise

def test_Kamodo_ambiguous_expr():
	a, b, c, x, y, z = symbols('a b c x y z')
	with pytest.raises(NameError):
		kamodo = Kamodo()
		kamodo[a(b,c)] = x**2 + y**2
	kamodo = Kamodo()
	kamodo['a(x, y)'] = x**2 + y**2
	assert kamodo.a(3,4) == 3**2 + 4**2

def test_Kamodo_callable():
	kamodo = Kamodo(f = lambda x: x**2, verbose = False)
	assert kamodo.f(3) == 9
	kamodo = Kamodo(f = lambda x,y: x**2 + y**3)
	assert kamodo.f(3,4) == 3**2 + 4**3


def test_Kamodo_callable_array():
	kamodo = Kamodo(f = lambda x: x**2)
	assert (kamodo.f(np.linspace(0,1,107)) == np.linspace(0,1,107)**2).all()


def test_Kamodo_str():
	kamodo = Kamodo('f(a,x,b) = a*x+b')
	try:
		assert kamodo.f(3,4,5) == 3*4+5
	except:
		print(kamodo.signatures)
		print(list(kamodo.keys()))
		raise

def test_Kamodo_latex():
	kamodo = Kamodo('$f(a,x,b) = a^x+b $')
	assert kamodo.f(3,4,5) == 3**4+5

def test_Kamodo_get():
	kamodo = Kamodo('$f(a,x,b) = a^x+b $')
	try:
		assert kamodo['f'](3,4,5) == 3**4+5
	except:
		print(kamodo.signatures)
		print(list(kamodo.keys()))
		raise


def test_Kamodo_set():
	kamodo = Kamodo()
	kamodo['f(a,x,b)'] = '$a^x+b$'
	assert kamodo.f(2,3,4) == 2**3+4


def test_Kamodo_mismatched_symbols():
	with pytest.raises(NameError):
		kamodo = Kamodo('$f(a,b) = a + b + c$', verbose = False)
		assert 'f' not in kamodo

def test_Kamodo_assign_by_expression():
	kamodo = Kamodo()
	f,x,y = symbols('f x y')
	kamodo[f(x,y)] = x**2 + y**2
	assert kamodo['f'](3,4) == 3**2 + 4**2
	assert kamodo.f(3,4) == 3**2 + 4**2 

def test_Kamodo_assign_by_callable():
	kamodo = Kamodo(f = lambda x: x**2, verbose = False)
	assert kamodo['f'](3) == 3**2
	assert kamodo.f(3) == 3**2	
	Kamodo(f = lambda x,y: x+y)
	def rho(x,y): 
		return x+y
	kamodo['g'] = rho
	assert kamodo.g(3,4) == 3 + 4

def test_Kamodo_composition():
	f, g = symbols('f g', cls = UndefinedFunction)
	x = Symbol('x')
	kamodo = Kamodo(f = "$x^2$")
	kamodo['g(x)'] = "$f(x) + x^2$"
	kamodo['h(x)'] = g(x)**2 + x**2

	assert kamodo.f(3) == 3**2
	assert kamodo.g(3) == 3**2 + 3**2
	assert kamodo.h(3) == (3**2 + 3**2)**2 + 3**2

def test_Kamodo_reassignment():
	a,b,c, x,y,z, r = symbols('a b c x y z r')
	kamodo = Kamodo('f(x) = 3*x+5')
	kamodo[r(x,y)] = x + y
	kamodo[r(x,y)] = "3*x + y"

def test_multivariate_composition():
	kamodo = Kamodo(f = 'x**2', g = lambda y: y**3)
	kamodo['h(x,y)'] = 'f(x) + g(y)'
	kamodo['i(x,y)'] = 'x + f(h(x,y))'
	assert kamodo.h(3,4) == 3**2 + 4**3
	assert kamodo.i(3,4) == 3 + (3**2 + 4**3)**2


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
	f = symbols('f', cls = UndefinedFunction)
	x = Symbol('x')
	f_ = f(x)
	kamodo = Kamodo()
	kamodo[f_] = x**2
	assert kamodo.f(3) == 9

def test_compact_variable_syntax():
	f = symbols('f', cls = UndefinedFunction)
	x = symbols('x') 
	kamodo = Kamodo(f = 'x**2')# f also registers f(x)
	kamodo['g(x)'] = 'f + x'
	assert kamodo.g(3) == 3**2 + 3

def test_unit_registry():
	def rho(x,y):
		return x*y
	def v(x,y):
		return x+y
	v.meta = dict(units = 'km/s')
	kamodo = Kamodo(v = v, verbose = False)
	kamodo['rho[kg/cc]'] = rho
	try:
		assert kamodo.rho.meta['units'] == 'kg/cc'
		assert kamodo.v.meta['units'] == 'km/s'
	except:
		print('\n',pd.DataFrame(kamodo.signatures))
		raise

	with pytest.raises(NameError):
		kamodo['p[kg]'] = v
	assert 'p' not in kamodo

def test_to_latex():
	kamodo = Kamodo(f = 'x**2')
	assert kamodo.to_latex() == '\\begin{equation}f{\\left (x \\right )} = x^{2}\\end{equation}'
	kamodo = Kamodo(g = 'x')
	assert kamodo.to_latex() == '\\begin{equation}g{\\left (x \\right )} = x\\end{equation}'


def test_expr_conversion():
	kamodo = Kamodo('$a[km] = x$')
	kamodo.a

def test_validate_units():
	f,x = symbols('f x')
	with pytest.raises(ValueError):
		validate_units(Eq(f,x), dict(f = get_unit('kg'), x = get_unit('m')))
	validate_units(Eq(f,x), dict(f = get_unit('kg'), x = get_unit('g')))
	
	lhs_units = validate_units(f(x), dict(f = get_unit('kg'), x = get_unit('m')))
	print(lhs_units)

def test_unit_conversion():
	kamodo = Kamodo('$a(x)[km/s] = x$',
						'$b(y)[m/s]  = y$',)
	kamodo['c(x,y)[m/s]'] = '$a + b$'
	assert kamodo.c(1,2) == 1002 


def test_get_unit():
	assert get_unit('kg/m^3') == get_unit('kg/m**3')

def test_unit_conversion_syntax():
	kamodo = Kamodo('rho[kg/m^3] = x', verbose = False)
	with pytest.raises(NameError):
		kamodo['d[kg]'] = 'rho'
		print(kamodo.detail())


def test_unit_composition():
	kamodo = Kamodo('m[kg] = x', verbose = False)
	kamodo['v[km/s]'] = 'y'
	kamodo['p'] = 'm*v'
	assert get_unit(kamodo.signatures['p(x, y)']['units']) == get_unit('kg*km/s')

def test_unit_function_composition():
	kamodo = Kamodo('X[m] = x', verbose = True)
	@kamodofy(units = 'km/s')
	def v(x):
		return y

	kamodo['v'] = v
	kamodo['speed'] = 'v(X)'
	assert kamodo.speed.meta['units'] == 'kilometer/second'



def test_method_args():
	class TestModel(Kamodo):
		def __init__(self, *args, **kwargs):
			super(TestModel, self).__init__(*args, **kwargs)
			self['rho'] = self.density
				
		def density(self, alt, lat, lon):
			return alt + lat + lon
		
		density.meta = dict(units = '1/cm^3')

	test = TestModel(verbose = False)
	assert str(list(test.keys())[0].args[0]) != 'self'

	test['r(alt, lat, lon)[1/m^3]'] = 'rho'
	try: # check that args are in the correct order
		assert list(test.signatures.values())[0]['symbol'].args == list(test.signatures.values())[1]['symbol'].args
	except:
		print('\n', test.detail())
		raise



def test_komodofy_decorator():
	@kamodofy(units = 'kg/cm^3')
	def my_density(x,y,z):
		return x+y+z

	assert my_density.meta['units'] == 'kg/cm^3'
	assert my_density(3,4,5) == 12

def test_komodofy_register():
	@kamodofy(units = 'kg/cm^3')
	def my_density(x,y,z):
		return x+y+z

	kamodo = Kamodo(rho = my_density)
	assert kamodo.rho.meta['units'] == 'kg/cm^3'

def test_komodofy_method():
	class TestModel(Kamodo):
		def __init__(self, *args, **kwargs):
			super(TestModel, self).__init__(*args, **kwargs)
			self['rho'] = self.density
		
		@kamodofy(units = '1/cm^3')        
		def density(self, alt, lat, lon):
			return alt + lat + lon

	test = TestModel(verbose = False)
	assert test.density.meta['units'] == '1/cm^3'



# 	try:
# 		assert len(kamodo.detail()) == 1
# 	except:
# 		print kamodo.symbol_registry
# 		print kamodo.detail()
# 		raise

# # def test_remove_symbol():
# # 	kamodo = Kamodo(f = 'x')
# # 	kamodo['g'] = '2*f'
# # 	kamodo.remove_symbol('f')
# # 	try:
# # 		assert len(kamodo.symbol_registry) == 1
# # 		assert len(kamodo) == 2
# # 	except:
# # 		print '\n',kamodo.detail()
# # 		raise

# def test_function_arg_ordering():
# 	def x(r,theta,phi):
# 	    '''converts from spherical to cartesian'''
# 	    return r*np.sin(theta)*np.cos(phi)

# 	kamodo = Kamodo(x = x)
# 	rhs_args = get_function_args(kamodo.x)
# 	lhs_args = kamodo.signatures.values()[0]['symbol'].args

# 	for i in range(3):
# 	    assert lhs_args[i] == rhs_args[i]

# def test_function_change_of_variables():
# 	def rho(x):
# 		return x**2

# 	kamodo = Kamodo(x = 'y')
# 	kamodo['rho'] = rho
# 	kamodo['rho_perp'] = 'rho(x)'

# 	assert str(get_function_args(kamodo.rho_perp)[0]) == 'y'


def test_kamodofy_update_attribute():
	@kamodofy(units = 'cm', update = 'y')
	def f(x):
		return x

	kamodo = Kamodo(f = f)
	assert f.update == 'y'


def test_kamodo_coupling():
	@kamodofy(units = 'cm', update = 'y')
	def y_iplus1(y, x):
		return y + x

	@kamodofy(units = 'm')
	def x_iplus1(x, y):
		return x - y


	kamodo = Kamodo()
	kamodo['x_iplus1'] = x_iplus1
	kamodo['y_iplus1'] = y_iplus1

	kamodo.x_iplus1.update = 'x'
	assert kamodo.x_iplus1.update == 'x'
	assert kamodo.y_iplus1.update == 'y'

	simulation = kamodo.simulate(y = 1, x = 0, steps = 1)

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
	@kamodofy(units = 'kg')
	def f(x, y):
		return x + y

	kamodo = Kamodo(f = f)
	kamodo.f([3,4,5],5)
	


def test_redefine_variable():
	kamodo = Kamodo(rho = 'x + y')
	kamodo['rho'] = 'a + b'
	kamodo['rho(a,b)'] = 'a*b'