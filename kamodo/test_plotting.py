from kamodo import Kamodo
import pytest
import numpy as np
import pandas as pd
from .plotting import scatter_plot, line_plot, vector_plot, contour_plot, surface, plane, tri_surface_plot, get_arg_shapes, plot_types, plot_dict, image, symbolic_shape


def test_scatter_plot():
    title = 'test plot title'
    kamodo_obj = Kamodo(f='x**3')
    result = kamodo_obj.evaluate('f', x=np.random.random((100, 3)))
    titles = dict(variable='x', title=title)

    [trace], plot_type, layout = scatter_plot(result, titles, True)

    assert trace['mode'] == 'markers'
    assert layout['title']['text'] == title
    assert plot_type == '3d-scatter'


def test_line_plot_line():
    kamodo = Kamodo(f='x**3')
    result = kamodo.evaluate('f', x=np.linspace(-1, 1, 100))

    titles = dict(variable='f',
                  title='$f(x)=x^3$',
                  title_lhs='f', title_short='f')

    traces, plot_type, layout = line_plot(result, titles, True)

    assert traces[0]['name'] == 'f'
    assert layout['title']['text'] == '$f(x)=x^3$'
    assert plot_type == 'line'


def test_line_plot_2d_line():
    def fvec(t):
        x = t * np.sin(t) ** 3
        y = t * np.cos(t) * 2
        return np.vstack((x, y)).T

    title = '$\\vec{f}(t)=sin(t)^3, cos(t)^2 \\; t \\in [-\\pi, \\pi]$'
    kamodo = Kamodo(fvec=fvec)
    result = kamodo.evaluate('fvec',
                             t=np.linspace(-np.pi, np.pi, 100))

    titles = dict(variable='fvec',
                  title=title,
                  title_lhs='fvec',
                  title_short='f',
                  units='')

    traces, plot_type, layout = line_plot(result, titles, True)

    assert plot_type == '2d-line'
    assert layout['title']['text'] == title


def test_line_plot_3d_line_pd():
    def fvec(t):
        x = t * np.sin(t) ** 3
        y = t * np.cos(t) * 2
        z = t * np.sin(t)
        return pd.DataFrame(dict(X_GSE=x, Y_GSE=y, Z_GSE=z))

    title = '$\\vec{f}(t)=sin(t)^3, cos(t)^2 \\; t \\in [-\\pi, \\pi]$'
    kamodo = Kamodo(fvec=fvec)
    results = kamodo.evaluate('fvec',
                              t=np.linspace(-np.pi, np.pi, 100))

    titles = dict(variable='fvec',
                  title=title,
                  title_lhs='fvec',
                  title_short='f',
                  units='')

    traces, plot_type, layout = line_plot(results, titles, True)

    assert plot_type == '3d-line'
    assert layout['title']['text'] == title


def test_vector_plot_2d_vector():
    x = np.linspace(-np.pi, np.pi, 25)
    y = np.linspace(-np.pi, np.pi, 30)
    xx, yy = np.meshgrid(x, y)
    points = np.array(list(zip(xx.ravel(), yy.ravel())))

    def fvec_Ncomma2(rvec_Ncomma2=points):
        ux = np.sin(rvec_Ncomma2[:, 0])
        uy = np.cos(rvec_Ncomma2[:, 1])
        return np.vstack((ux, uy)).T

    title = '$\\vec{f}_{N,2}{\\left(\\vec{r}_{N,2} \\right)} = \\lambda{\\left(\\vec{r}_{N,2} \\right)}$'
    kamodo = Kamodo(fvec_Ncomma2=fvec_Ncomma2)
    results = kamodo.evaluate('fvec_Ncomma2', rvec_Ncomma2=points)

    titles = dict(variable='fvec_Ncomma2', title=title)
    [trace], plot_type, layout = vector_plot(results, titles, True, scale=0.1)

    assert trace['mode'] == 'lines'
    assert trace['name'] == 'quiver'
    assert plot_type == '2d-vector'
    assert layout['title']['text'] == title


def test_vector_plot_3d_vector():
    x, y, z = np.meshgrid(np.linspace(-2, 2, 4),
                          np.linspace(-3, 3, 6),
                          np.linspace(-5, 5, 10))
    points = np.array(list(zip(x.ravel(), y.ravel(), z.ravel())))

    def fvec_Ncomma3(rvec_Ncomma3=points):
        return rvec_Ncomma3

    title = '$\\vec{f}_{N,3}{\\left(\\vec{r}_{N,3} \\right)} = \\lambda{\\left(\\vec{r}_{N,3} \\right)}$'
    kamodo = Kamodo(fvec_Ncomma3=fvec_Ncomma3)

    results = kamodo.evaluate('fvec_Ncomma3', rvec_Ncomma3=points)

    titles = dict(variable='fvec_Ncomma3', title=title)
    [trace], plot_type, layout = vector_plot(results, titles, True)

    assert trace['hoverinfo'] == 'x+y+z+u+v+w+norm'
    assert plot_type == '3d-vector'
    assert layout['title']['text'] == title


def test_vector_plot_3d_vector_pd():
    x, y, z = np.meshgrid(np.linspace(-2, 2, 4),
                          np.linspace(-3, 3, 6),
                          np.linspace(-5, 5, 10))
    points = pd.DataFrame(list(zip(x.ravel(), y.ravel(), z.ravel())))

    def fvec_Ncomma3(rvec_Ncomma3=points):
        return rvec_Ncomma3

    title = '$\\vec{f}_{N,3}{\\left(\\vec{r}_{N,3} \\right)} = \\lambda{\\left(\\vec{r}_{N,3} \\right)}$'
    kamodo = Kamodo(fvec_Ncomma3=fvec_Ncomma3)

    results = kamodo.evaluate('fvec_Ncomma3', rvec_Ncomma3=points)

    titles = dict(variable='fvec_Ncomma3', title=title)

    [trace], plot_type, layout = vector_plot(results, titles, True)

    assert trace['hoverinfo'] == 'x+y+z+u+v+w+norm'
    assert plot_type == '3d-vector'
    assert layout['title']['text'] == title


def test_contour_plot_2d_grid():
    def f_NcommaM(x_N=np.linspace(0, 8 * np.pi, 100), y_M=np.linspace(0, 5, 90)):
        x, y = np.meshgrid(x_N, y_M, indexing='xy')
        return np.sin(x) * y

    title = '$\\operatorname{f_{N,M}}{\\left(x_{N},y_{M} \\right)} = \\lambda{\\left(x_{N},y_{M} \\right)}$'
    kamodo = Kamodo(f_NcommaM=f_NcommaM)

    results = kamodo.evaluate('f_NcommaM')
    titles = dict(variable='f_NcommaM', title=title)

    traces, chart_type, layout = contour_plot(results, titles, 'ij', True)

    assert chart_type == '2d-grid'
    assert layout['title']['text'] == title


def test_coutour_plot_ij_index_2d_grid():
    def f_NcommaM(x_N=np.linspace(0, 8 * np.pi, 100), y_M=np.linspace(0, 5, 100)):
        x, y = np.meshgrid(x_N, y_M, indexing='xy')
        return np.sin(x) * y

    title = '$\\operatorname{f_{N,M}}{\\left(x_{N},y_{M} \\right)} = \\lambda{\\left(x_{N},y_{M} \\right)}$'
    kamodo = Kamodo(f_NcommaM=f_NcommaM)

    results = kamodo.evaluate('f_NcommaM')
    titles = dict(variable='f_NcommaM', title=title)

    traces, chart_type, layout = contour_plot(results, titles, 'ij', True)

    assert chart_type == '2d-grid'
    assert layout['title']['text'] == title


def test_coutour_plot_xy_index_2d_grid():
    def f_NcommaM(x_N=np.linspace(0, 8 * np.pi, 100), y_M=np.linspace(0, 5, 100)):
        x, y = np.meshgrid(x_N, y_M, indexing='xy')
        return np.sin(x) * y

    title = '$\\operatorname{f_{N,M}}{\\left(x_{N},y_{M} \\right)} = \\lambda{\\left(x_{N},y_{M} \\right)}$'
    kamodo = Kamodo(f_NcommaM=f_NcommaM)

    results = kamodo.evaluate('f_NcommaM')
    titles = dict(variable='f_NcommaM', title=title)

    traces, chart_type, layout = contour_plot(results, titles, 'xy', True)

    assert chart_type == '2d-grid'
    assert layout['title']['text'] == title


def test_contour_plot_2d_skew():
    r = np.linspace(1, 3, 20)
    theta = np.linspace(0, np.pi, 14)
    r_, theta_ = np.meshgrid(r, theta)
    XX = r_ * np.cos(theta_)
    YY = r_ * np.sin(theta_)

    def f_NM(x_NM=XX, y_NM=YY):
        return np.sin(x_NM) + y_NM

    title = '$\\operatorname{f_{NM}}{\\left(x_{NM},y_{NM} \\right)} = \\lambda{\\left(x_{NM},y_{NM} \\right)}$'
    kamodo = Kamodo(f_NM=f_NM)

    results = kamodo.evaluate('f_NM')

    titles = dict(variable='f_NM', title=title)

    traces, chart_type, layout = contour_plot(results, titles, 'xy', True)

    assert chart_type == '2d-skew'
    assert layout['title']['text'] == title


def test_plane():
    def f_LMN(
            x_L=np.linspace(-5, 5, 50),
            y_M=np.linspace(0, 10, 75),
            z_N=np.linspace(-20, 20, 100)):
        xx, yy, zz = np.meshgrid(x_L, y_M, z_N, indexing='xy')
        return xx + yy + zz

    title = '$\\operatorname{f_{LMN}}{\\left(x_{L},y_{M},z_{N} \\right)} = \\lambda{\\left(x_{L},y_{M},z_{N} \\right)}$'
    kamodo = Kamodo(f_LMN=f_LMN)
    results = kamodo.evaluate('f_LMN', x_L=-5)
    titles = dict(variable='f_LMN', title=title)
    traces, chart_type, layout = plane(results, titles, verbose=True)

    assert chart_type == '3d-plane'
    assert layout['title']['text'] == title


def test_plane_meshgrid_false():
    def f_LMN(
            x_L=pd.DataFrame(np.linspace(-5, 5, 50)),
            y_M=pd.DataFrame(np.linspace(0, 10, 75)),
            z_N=pd.DataFrame(np.linspace(-20, 20, 75))):
        xx, yy, zz = np.meshgrid(x_L, y_M, z_N, indexing='xy')
        return xx + yy + zz

    title = '$\\operatorname{f_{LMN}}{\\left(x_{L},y_{M},z_{N} \\right)} = \\lambda{\\left(x_{L},y_{M},z_{N} \\right)}$'
    kamodo = Kamodo(f_LMN=f_LMN)
    results = kamodo.evaluate('f_LMN', x_L=-5)
    titles = dict(variable='f_LMN', title=title)

    traces, chart_type, layout = plane(results, titles, verbose=True)

    assert chart_type == '3d-plane'
    assert layout['title']['text'] == title


def test_surface_3d_surface():
    u = np.linspace(-2, 2, 40)
    v = np.linspace(-2, 2, 50)
    uu, vv = np.meshgrid(u, v)

    def parametric(x_NM=uu * np.sin(vv * np.pi),
                   y_NM=vv,
                   z_NM=np.exp(-uu ** 2 - vv ** 2)):
        return np.array([1])

    title = '$p{\\left(x_{NM},y_{NM},z_{NM} \\right)} = \\lambda{\\left(x_{NM},y_{NM},z_{NM} \\right)}$'
    kamodo = Kamodo(p=parametric)

    results = kamodo.evaluate('p')
    titles = dict(variable='p', title=title)

    traces, chart_type, layout = surface(results, titles)

    assert chart_type == '3d-surface'
    assert layout['title']['text'] == title


def test_surface_3d_surface_scalar():
    R = 1
    theta = np.linspace(.2 * np.pi, .8 * np.pi, 40)
    phi = np.linspace(0, 2 * np.pi, 50)
    theta_, phi_ = np.meshgrid(theta, phi)
    r = (R + .1 * (np.cos(10 * theta_) * np.sin(14 * phi_)))

    xx = r * np.sin(theta_) * np.cos(phi_)
    yy = r * np.sin(theta_) * np.sin(phi_)
    zz = r * np.cos(theta_)

    def spherelike(x_NM=xx, y_NM=yy, z_NM=zz):
        return .1 * x_NM + x_NM ** 2 + y_NM ** 2 + z_NM ** 2

    title = '$\\operatorname{h_{NM}}{\\left(x_{NM},y_{NM},z_{NM} \\right)} = \\lambda{\\left(x_{NM},y_{NM},z_{NM} \\right)}$'
    kamodo = Kamodo(h_NM=spherelike)

    results = kamodo.evaluate('h_NM')
    titles = dict(variable='h_NM', title=title)

    traces, chart_type, layout = surface(results, titles, True)

    assert chart_type == '3d-surface-scalar'
    assert layout['title']['text'] == title


def test_surface_len_of_result_variable_is_1():
    def f_LMN(
            x_L=np.linspace(-5, 5, 50),
            y_M=np.linspace(0, 10, 1),
            z_N=np.linspace(-20, 20, 1)):
        xx, yy, zz = np.meshgrid(x_L, y_M, z_N, indexing='xy')
        return xx + yy + zz

    title = '$\\operatorname{f_{LMN}}{\\left(x_{L},y_{M},z_{N} \\right)} = \\lambda{\\left(x_{L},y_{M},z_{N} \\right)}$'
    kamodo = Kamodo(f_LMN=f_LMN)
    results = kamodo.evaluate('f_LMN')
    titles = dict(variable='f_LMN', title=title)

    with pytest.raises(UnboundLocalError) as error:
        traces, chart_type, layout = surface(results, titles, True)

    assert "local variable 'chart_type' referenced before assignment" in str(error)


def test_tri_surface_plot():
    R = 1
    theta = np.linspace(.2 * np.pi, .8 * np.pi, 40)
    phi = np.linspace(0, 2 * np.pi, 50)
    theta_, phi_ = np.meshgrid(theta, phi)
    r = (R + .1 * (np.cos(10 * theta_) * np.sin(14 * phi_)))

    xx = r * np.sin(theta_) * np.cos(phi_)
    yy = r * np.sin(theta_) * np.sin(phi_)
    zz = r * np.cos(theta_)

    def tri_surface(x_NM=xx, y_NM=yy, z_NM=zz):
        return .1 * x_NM + x_NM ** 2 + y_NM ** 2 + z_NM ** 2

    title = '$\\operatorname{h_{NM}}{\\left(x_{NM},y_{NM},z_{NM} \\right)} = \\lambda{\\left(x_{NM},y_{NM},' \
            'z_{NM} \\right)}$ '
    kamodo = Kamodo(h_NM=tri_surface)

    results = kamodo.evaluate('h_NM')
    titles = dict(variable='h_NM', title=title)

    traces, chart_type, layout = tri_surface_plot(results, titles, True)

    assert chart_type == '3d-surface'
    assert layout['title']['text'] == title


def test_arg_shape_np():
    def fvec(t):
        x = t * np.sin(t) ** 3
        y = t * np.cos(t) * 2
        return np.vstack((x, y)).T

    kamodo = Kamodo(fvec=fvec)
    results = kamodo.evaluate('fvec',
                              t=np.linspace(-np.pi, np.pi, 100))

    f = results['fvec']

    shape = get_arg_shapes(f)

    assert type(f) == np.ndarray
    assert shape == [(100, 2)]


def test_arg_shape_pd():
    def fvec(t):
        x = t * np.sin(t) ** 3
        y = t * np.cos(t) * 2
        return pd.DataFrame(dict(X_GSE=x, Y_GSE=y))

    kamodo = Kamodo(fvec=fvec)
    results = kamodo.evaluate('fvec',
                              t=np.linspace(-np.pi, np.pi, 100))

    f = results['fvec']

    shape = get_arg_shapes(f)

    assert type(f) == pd.DataFrame
    assert shape == [(100, 2)]


def test_line_plot_len_of_result_3():
    r = np.linspace(1, 3, 20)
    theta = np.linspace(0, np.pi, 14)
    r_, theta_ = np.meshgrid(r, theta)
    XX = r_ * np.cos(theta_)
    YY = r_ * np.sin(theta_)

    def f_NM(x_NM=XX, y_NM=YY):
        return np.sin(x_NM) + y_NM

    title = '$\\operatorname{f_{NM}}{\\left(x_{NM},y_{NM} \\right)} = \\lambda{\\left(x_{NM},y_{NM} \\right)}$'
    kamodo = Kamodo(f_NM=f_NM)
    titles = dict(variable='f_NM', title=title)

    results = kamodo.evaluate('f_NM')

    [trace], chart_type, layout = line_plot(results, titles, True)

    assert trace['mode'] == 'lines'
    assert chart_type == '2d-line'
    assert layout['title'] == title


def test_line_plot_len_of_result_4():
    r = np.linspace(1, 3, 20)
    theta = np.linspace(0, np.pi, 14)
    r_, theta_ = np.meshgrid(r, theta)
    XX = r_ * np.cos(theta_)
    YY = r_ * np.sin(theta_)
    ZZ = r_ * np.sin(theta_)

    def f_NM(x_NM=XX, y_NM=YY, z_NM=ZZ):
        return np.sin(x_NM) + y_NM + z_NM

    title = '$\\operatorname{f_{NM}}{\\left(x_{NM},y_{NM},z_{NM} \\right)} = \\lambda{\\left(x_{NM},y_{NM},z_{NM} \\right)}$'
    kamodo = Kamodo(f_NM=f_NM)

    titles = dict(variable='f_NM', title=title)
    results = kamodo.evaluate('f_NM')

    [trace], chart_type, layout = line_plot(results, titles, True)

    assert trace['mode'] == 'lines'
    assert chart_type == '4d-line'
    assert layout['title']['text'] == title


def test_scatter_plot_pd():
    kamodo = Kamodo(f='x**3')
    r = pd.DataFrame({'x': np.linspace(0, 1, 100), 'y': np.linspace(-1, 1, 100), 'z': np.linspace(-1, 0, 100)})

    results = kamodo.evaluate('f', x=r)
    title = '$f{\\left(x \\right)} = x^{3}$'
    titles = dict(variable='f',
                  title=title,
                  title_lhs='f', title_short='f')

    [trace], plot_type, layout = scatter_plot(results, titles, True)

    assert trace['mode'] == 'markers'
    assert plot_type == '3d-scatter'
    assert layout['title']['text'] == title


def test_line_plot_3d_line():
    def fvec(t):
        x = t * np.sin(t) ** 3
        y = t * np.cos(t) * 2
        z = t * np.cos(t)
        return np.vstack((x, y, z)).T

    title = '$\\vec{f}{\\left(t \\right)} = \\lambda{\\left(t \\right)}$'
    kamodo = Kamodo(fvec=fvec)
    results = kamodo.evaluate('fvec',
                              t=np.linspace(-np.pi, np.pi, 100))

    titles = dict(variable='fvec',
                  title=title,
                  title_lhs='fvec',
                  title_short='f',
                  units='')

    [trace], chart_type, layout = line_plot(results, titles, True)
    assert trace['mode'] == 'lines'
    assert chart_type == '3d-line'
    assert layout['title']['text'] == title

def test_image_plot():
    def img(i=np.arange(33),j=np.arange(35)):
        ii, jj, kk = np.meshgrid(i,j,[100, 200, 255], indexing='ij')
        return kk

    kamodo = Kamodo(img=img)
    results = kamodo.evaluate('img')
    titles = dict(variable='img', title='mytitle')
    [trace], chart_type, layout = image(results, titles, True)


def test_plot_keys():
    for k in plot_types.to_dict(orient = 'index'):
        try:
            plot_dict[k[0]][k[1]]['name']
        except KeyError:
            print('could not find', k[0], k[1])
            raise


def test_symbolic_shape():
    assert symbolic_shape((3,1,3)) == ((3, 1, 3),)
    assert symbolic_shape((4,4,4)) == (('N', 'N', 'N'),)
    assert symbolic_shape((2,1,3)) == ((2, 1, 3),)
    assert symbolic_shape((4,5,6)) == (('N', 'M', 'L'),)
    assert symbolic_shape((4,5,1)) == (('N', 'M', 1),)
    assert symbolic_shape((5, 6), (6, 5), (5, 4)) == (('N', 'M'), ('M', 'N'), ('N', 'L'))

