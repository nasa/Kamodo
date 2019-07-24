"""
Copyright © 2017 United States Government as represented by the Administrator, National Aeronautics and Space Administration.  
No Copyright is claimed in the United States under Title 17, U.S. Code.  All Other Rights Reserved.
"""
import plotly.graph_objs as go
import numpy as np
from util import arg_to_latex, beautify_latex, cast_0_dim
from plotly import figure_factory as ff
import pandas as pd
from collections import defaultdict


def line_plot(result, titles, verbose = False, **kwargs):
	'''N-d line plot f(t)'''
	f = result[titles['variable']]
	t_name, t = list(result.items())[0]
	if len(result) == 2:
		if len(f.shape) == 1:
			if verbose:
				print('\t1-d output', f.shape)
			assert f.shape == t.shape
			y = f
			layout = go.Layout(
				title = titles['title'], 
				xaxis = dict(title = '${}$'.format(t_name)),
				yaxis = dict(title = titles['title_lhs']))
			trace = go.Scatter(x = t, y = y, name = titles['title_short'])
			chart_type = 'line'

		elif f.shape[1] == 2:
			if verbose:
				print('\t2-d output', f.shape)
			x = f[:,0]
			y = f[:,1]
			# t_name, t = result.items()[0]
			text = ["{}:{}".format(t_name, v) for v in t]
			layout = go.Layout(title = titles['title'], 
							   xaxis = dict(title = 'x' + titles['units']),
							  yaxis = dict(title = 'y' + titles['units']))
			trace = go.Scatter(x = x, y = y, text = text)
			chart_type = '2d-line'

		elif f.shape[1] == 3: 
			if verbose:
				print('\t3-d output', f.shape)
			x = f[:,0]
			y = f[:,1]
			z = f[:,2]
			# t_name, t = result.items()[0]
			text = ["{}:{}".format(t_name, v) for v in t]
			layout = go.Layout(
				title = titles['title'],
				scene = dict(
					xaxis = dict(title = 'x' + titles['units']),
					yaxis = dict(title = 'y' + titles['units']),
					zaxis = dict(title = 'z' + titles['units']),
					))
			trace = go.Scatter3d(x = x, y = y, z = z, text = text, mode = 'lines')
			chart_type = '3d-line'

	elif len(result) == 4:
		arg0, val0 = t_name, t
		arg1, val1 = list(result.items())[1]
		arg2, val2 = list(result.items())[2]

		assert(f.shape == t.shape)

		text = ["{}:{}".format(titles['variable'], v) for v in f]
		
		trace = go.Scatter3d(
			x = val0, 
			y = val1, 
			z = val2, 
			text = text,
			marker = dict(
				color = result[titles['variable']], 
				colorscale = 'Viridis',
				colorbar = dict(title = titles['title_short'])),
			line = dict(
				color = result[titles['variable']],
				colorscale='Viridis'), 
			mode = 'markers+lines')
		layout = go.Layout(
			title = titles['title'],
			scene = dict(
				xaxis = dict(title = arg0),
				yaxis = dict(title = arg1),
				zaxis = dict(title = arg2),))
		chart_type = '4d-line'
	else:
		raise NotImplementedError('shape not supported: {}'.format(f.shape))

	return [trace], chart_type, layout


def vector_plot(result, titles, verbose = False, **kwargs):
	variable = titles['variable']
	val0 = list(result.values())[0]

	if (result[variable].shape == val0.shape) & (val0.shape[1] == 2):
		if verbose:
			print('\t 2-d output', result[variable].shape)
		u = result[variable][:,0]
		v = result[variable][:,1]
		x = val0[:,0]
		y = val0[:,1]

		trace = ff.create_quiver(
			x, y, u, v,
			# scale=.25,
			# arrow_scale=.4,
			name='quiver',
			line=dict(width=1))['data'][0] # extracts quiver trace

		layout = go.Layout(title = titles['title'], 
							  xaxis = dict(title = 'x'),
							  yaxis = dict(title = 'y'))
		chart_type = '2d-vector'

	elif (result[variable].shape == val0.shape) & (val0.shape[1] == 3):
		if verbose:
			print('\t 3-d output', result[variable].shape)
		u = result[variable][:,0].tolist()
		v = result[variable][:,1].tolist()
		w = result[variable][:,2].tolist()
		x = val0[:,0].tolist()
		y = val0[:,1].tolist()
		z = val0[:,2].tolist()
		trace = go.Cone(x = x, y = y, z = z, u = u, v = v, w = w)

		layout = go.Layout(
			title = titles['title'],
			scene = dict(
				xaxis = dict(title = 'x'),
				yaxis = dict(title = 'y'),
				zaxis = dict(title = 'z'),
				))
		chart_type = '3d-vector'
	else:
		raise NotImplementedError('plot type not supported yet')

	return [trace], chart_type, layout

def contour_plot(result, titles, indexing, verbose = False, **kwargs):
	variable = titles['variable']
	if verbose:
		print('\t2-d output', result[variable].shape)
	z = result[variable]
	title = titles['title']
	arg0, val0 = list(result.items())[0]
	arg1, val1 = list(result.items())[1]

	if verbose:
		print('2 inputs: {} {}, {} {}'.format(arg0, val0.shape, arg1, val1.shape))
	traces = []
	if len(val0.shape) == 1:
		if verbose:
			print('\t\t1-d args', end=' ')
		if result[variable].shape[0] == val0.shape[0]:
			if indexing == 'ij':
				if verbose:
					print('{} indexing'.format(indexing))
				trace = go.Contour(x = val0,
								   y = val1,
								   z = z.T)
			else:
				if verbose:
					print('{} indexing'.format(indexing))
				trace = go.Contour(x = val0,
								   y = val1,
								   z = z)

		else:
			if verbose:
				print('xy indexing')
			trace = go.Contour(x = val0,
							   y = val1,
							   z = z)
		layout = go.Layout(
			title = title,
			xaxis = dict(title = '${}$'.format(arg0)),
			yaxis = dict(title = '${}$'.format(arg1)))
		traces.append(trace)
		chart_type = '2d-grid'
	else:
		if verbose:
			print('\t\t2-d args', val0.shape, val1.shape)
		assert val0.shape == val1.shape
		assert val0.shape == z.shape
		xaxis = dict(title = beautify_latex('${}$'.format(arg_to_latex(arg0))))
		yaxis = dict(title = beautify_latex('${}$'.format(arg_to_latex(arg1))))
		carpet_traces, layout = carpet_plot(result, title, xaxis, yaxis)
		traces.extend(carpet_traces)
		chart_type = '2d-skew'

	return traces, chart_type, layout

def carpet_plot(results, title, xaxis, yaxis, indexing = 'xy', **kwargs):
	'''Assumes ordered dict where values have the same shape'''
	# print title, xaxis, yaxis
	arg0, val0 = list(results.items())[0]
	arg1, val1 = list(results.items())[1]
	arg2, val2 = list(results.items())[2]
	
	a = list(range(val0.shape[0]))
	b = list(range(val0.shape[1]))

	aa, bb = np.meshgrid(b,a, indexing = indexing)

	trace1 = go.Contourcarpet(
		a = aa.ravel(),
		b = bb.ravel(),
		z = val2.ravel(),
	)

	trace2 = go.Carpet(
		a = aa.ravel(),
		b = bb.ravel(),
		x = val0.ravel(),
		y = val1.ravel(),
	)
	traces = [trace1, trace2]
	layout = go.Layout(
		title = title, 
		xaxis = xaxis,
		yaxis =  yaxis)
	return traces, layout

def plane(result, titles, indexing = 'xy', verbose = False, **kwargs):
	variable = titles['variable']
	arg0, val0 = list(result.items())[0]
	arg1, val1 = list(result.items())[1]
	arg2, val2 = list(result.items())[2]
	traces = []
	if verbose:
		print('\t3-d output', result[variable].shape)
	surfacecolor = result[variable].squeeze()
	create_meshgrid = True
	for x in val0, val1, val2:
		if len(get_arg_shapes(x)[0]) != 1:
			create_meshgrid = False
			ones = np.ones_like(x)
	if create_meshgrid:
		xx, yy, zz = np.meshgrid(val0, val1, val2)
	else:
		xx, yy, zz = [v*ones for v in [val0, val1, val2]]
	traces.append(go.Surface(
		x = xx.squeeze(),
		y = yy.squeeze(),
		z = zz.squeeze(),
		surfacecolor = surfacecolor))
	layout = go.Layout(
		title = titles['title'],
		scene = dict(
			xaxis = dict(title = arg0),
			yaxis = dict(title = arg1),
			zaxis = dict(title = arg2),
		)
	)
	chart_type = '3d-plane'

	return traces, chart_type, layout

def surface(result, titles, verbose = False, **kwargs):
	variable = titles['variable']
	title = titles['title']
	arg0, val0 = list(result.items())[0]
	arg1, val1 = list(result.items())[1]
	arg2, val2 = list(result.items())[2]

	if verbose:
		print('\t2-d output', result[variable].shape)

	traces = []

	if len(result[variable].shape) == 2:
		surfacecolor = result[variable]
		xx = cast_0_dim(val0, result[variable])
		yy = cast_0_dim(val1, result[variable])
		zz = cast_0_dim(val2, result[variable])
		traces.append(go.Surface(
			x = xx,
			y = yy,
			z = zz,
			surfacecolor = surfacecolor))
		layout = go.Layout(
			title = title,
			scene = dict(
				xaxis = dict(title = arg0),
				yaxis = dict(title = arg1),
				zaxis = dict(title = arg2),))
		chart_type = '3d-surface-scalar'

	elif len(result[variable].shape) == 1:
		if verbose:
			print('\tscalar output', result[variable].shape)

		elif result[variable].shape[0] == 1:
			traces.append(go.Surface(
				x = val0, 
				y = val1, 
				z = val2, 
				))
			layout = go.Layout(
				title = title,
				scene = dict(
					xaxis = dict(title = arg0),
					yaxis = dict(title = arg1),
					zaxis = dict(title = arg2),
					))
			chart_type = '3d-surface'

	elif len(result[variable]) == 1: #scalar or color value?
		if verbose:
			print('\tsingle valued output', result[variable].shape)
		surfacecolor = result[variable]
		arg0, val0 = list(result.items())[0]
		arg1, val1 = list(result.items())[1]
		arg2, val2 = list(result.items())[2]
		xx = val0
		yy = val1
		zz = val2
		traces.append(go.Surface(
			x = xx,
			y = yy,
			z = zz,
			))
		layout = go.Layout(
			title = title,
			scene = dict(
				xaxis = dict(title = arg0),
				yaxis = dict(title = arg1),
				zaxis = dict(title = arg2),
			)
		)

	return traces, chart_type, layout



plot_dict = {
	(1,)	:	{(('N','M'), ('N','M'), ('N','M')): {'name': '3d-parametric', 'func': surface}},
	('N',)	:	{
		(('N',),) : {'name': '1d-line', 'func': line_plot},
		(('N',),('N',),('N',)): {'name': '3d-line-scalar', 'func': line_plot},},
	('N',2) :	{
		(('N',),) :{'name': '2d-line', 'func': line_plot},
		(('N',2),):{'name': '2d-vector', 'func': vector_plot},
	},
	('N',3) : 	{
		(('N',),) :{'name': '3d-line', 'func': line_plot},
		(('N',3),):{'name': '3d-vector', 'func': vector_plot},
	},
	('N','M'):	{
		(('N',),('M',)) :{'name': '2d-contour', 'func': contour_plot},
		(('N','M'),('N','M')):{'name': '2d-contour-skew', 'func': contour_plot},
		(('N','M'), ('N','M'), ('N','M')): {'name': '3d-parametric-scalar', 'func': surface},
		((1,),('N','M'),('N','M')):{'name': '3d-plane', 'func': plane},
		(('N','M'),(1,),('N','M')):{'name': '3d-plane', 'func': plane},
		(('N','M'),('N','M'),(1,)):{'name': '3d-plane', 'func': plane},
	},
	('N','M',1): {
		((1,),('N',),('M',)):{'name': '3d-plane', 'func': plane},
		(('N',),(1,),('M',)):{'name': '3d-plane', 'func': plane},
		(('N',),('M',),(1,)):{'name': '3d-plane', 'func': plane},
	},
}

def get_arg_shapes(*args):
	shapes = []
	for a in args:
		if type(a) == np.ndarray:
			shape = a.shape
		else:
			try:
				shape = (len(a),)
			except:
				shape = (1,)
		shapes.append(shape)
	return shapes


def get_plot_key(out_shape, *arg_shapes):
	shapes_match = all([(a == out_shape) for a in list(arg_shapes)])
	nargs = len(arg_shapes)
	arg_dims = ''
	if nargs == 1:
		if len(out_shape) == 1:
			out_dim = 'N'
		else:
			out_dim = ('N', out_shape[-1])
	else:
		if len(out_shape) == 1:
			if out_shape[0] == 1:
				out_dim = 1,
			else:
				out_dim = 'N',
		elif len(out_shape) == 2:
			out_dim = 'N','M'
		elif len(out_shape) == 3:
			if 1 in out_shape:
				out_dim = 'N', 'M', 1
			else:
				out_dim = 'N','M','L'
	
	out_dim = tuple(out_dim)
	if shapes_match:
		arg_dims = tuple(nargs*[out_dim])
	else:
		if nargs == 1:
			if out_dim == ('N',2):
				if arg_shapes[0][0] == out_shape[0]:
					arg_dims = (('N',),)
			elif out_dim == ('N',3):
				if arg_shapes[0][0] == out_shape[0]:
					arg_dims = (('N',),)
		elif nargs == 2:
			if out_dim == ('N','M'):
				if (arg_shapes[0][0] in out_shape) & (arg_shapes[1][0] in out_shape):
					arg_dims = (('N',),('M',))
		elif nargs == 3:
			if out_dim == (1,):
				if len(set(arg_shapes)) == 1:
					if len(arg_shapes[0]) == 2:
						arg_dims = tuple(3*[('N','M')])
			elif out_dim == ('N','M'):
				arg_set = set([1])
				for arg in arg_shapes:
					arg_set.update(set(arg))
				if arg_set - set(out_shape) == set([1]):
					if arg_shapes[0] == (1,):
						arg_dims = ((1,), ('N','M'), ('N','M'))
					elif arg_shapes[1] == (1,):
						arg_dims = (('N','M'),(1,),('N','M',))
					elif arg_shapes[2] == (1,):
						arg_dims = (('N','M'),('N','M'),(1,))
			elif out_dim == ('N','M',1):
				arg_set = set([1])
				for arg in arg_shapes:
					arg_set.update(set(arg))
				if arg_shapes[0] == (1,):
					arg_dims = ((1,), ('N',), ('M',))
				elif arg_shapes[1] == (1,):
					arg_dims = (('N',),(1,),('M',))
				elif arg_shapes[2] == (1,):
					arg_dims = (('N',),('M',),(1,))
	if arg_dims == '':
		raise NotImplementedError('No way to handle out_shape {}, with arg shapes:{}'.format(out_shape, arg_shapes))
	return out_dim, arg_dims


plot_types = dict()
for out_shape, v in list(plot_dict.items()):
	for arg_shapes, v_ in list(v.items()):
		plot_key = get_plot_key(out_shape,*arg_shapes)
		if plot_key in plot_types:
			raise KeyError('plot_key already present {}'.format(plot_key))
		plot_types[plot_key] = [v_['name'], v_['func']]

plot_types = pd.DataFrame(plot_types).T
plot_types.index.set_names(['out_shape', 'arg_shapes'], inplace = True)
plot_types.columns = ['plot_type', 'function']

def test_plot_keys():
	for k in plot_types.to_dict(orient = 'index'):
		try:
			plot_dict[k[0]][k[1]]['name']
		except KeyError:
			print('could not find', k[0], k[1])
			raise

