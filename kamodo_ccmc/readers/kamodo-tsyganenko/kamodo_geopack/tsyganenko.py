# +
from kamodo import Kamodo, kamodofy
import PyGeopack as gp
import numpy as np
import pandas as pd
import scipy

from datetime import datetime, timezone
from kamodo.plotting import plot_types


# +
def to_utc(t_):
    """inputs is assumed to be a list of times"""
    t = pd.Series(t_)
    if t.dt.tz is None:
        t = t.dt.tz_localize('UTC')
    else:
        t = t.dt.tz_convert('UTC')
    return t

def time_to_tsyg(time):
    """time must be pandas time series
    returns two numpy arrays for date and time"""
    date_int = time.dt.strftime('%Y%m%d').astype(int)
    ut = time.dt.hour + time.dt.minute/60 + (time.dt.second + 1e-6*time.dt.microsecond)/3600
    if len(time) > 1:
        return date_int.to_numpy(), ut.to_numpy()
    else:
        return date_int.to_numpy()[0], ut.to_numpy()[0]

def get_component(obj, component, coord):
    return getattr(obj, component + coord)
    

class KTsyganenko(Kamodo): 
    
    def get_model_params(self, t):
        time = to_utc(t)
        date, ut = time_to_tsyg(time)
        return gp.Params.GetModelParams(date, ut, Model = self.model_name)

    def trace(self, x, y, z, date_int, ut,
                 model_name=None,
                 coord_in=None,
                 coord_out=None):
        """pass through for gp.ModelField"""
        if model_name is None:
            model_name = self.model_name
        if coord_in is None:
            coord_in = self.coord_in
        if coord_out is None:
            coord_out = self.coord_out

        trace = gp.TraceField(x, y, z, date_int, ut,
                                 Model=model_name,
                                 CoordIn=coord_in,
                                 CoordOut=coord_out)
        return trace
    
    def __init__(self, model_name='TS05', coord_in='GSM', coord_out='GSM', **kwargs):
        
        self.model_name = model_name
        self.coord_in = coord_in
        self.coord_out = coord_out

        # Prepare model for function registration for the input argument
        super(KTsyganenko, self).__init__(**kwargs)


        in_vec_sym = '\\vec{x}_{' + self.coord_in + '}'
        out_vec_sym = '\\vec{B}_{' + self.coord_out + '}'
        x_in_sym = 'x_{' + self.coord_in + '}'
        y_in_sym = 'y_{' + self.coord_in + '}'
        z_in_sym = 'z_{' + self.coord_in + '}'

        @kamodofy(
            units = 'nT',
            equation=f'{out_vec_sym}({x_in_sym}, {y_in_sym}, {z_in_sym}, t)',
            citation=f'Tsynganenko version: {model_name}. \nFor model details: https://ccmc.gsfc.nasa.gov/models/Tsyganenko%20Magnetic%20Field~TS05/ \nBuilt on PyGeopack https://github.com/mattkjames7/PyGeopack')
        def Bfield(x, y, z, t):
            """Magnetic field model from Tsyganenko
            see self.model_name for specific model used
            """
            time = to_utc(t)
            date_int, ut = time_to_tsyg(time)
            b = gp.ModelField(x, y, z, date_int, ut,
                                 Model=self.model_name,
                                 CoordIn=self.coord_in,
                                 CoordOut=self.coord_out)
            return np.hstack(b)
        

        @kamodofy(
            units = 'nT',
            equation= f'{out_vec_sym}({in_vec_sym}, t)',
            citation=f'Tsynganenko version: {model_name}. \nFor model details: https://ccmc.gsfc.nasa.gov/models/Tsyganenko%20Magnetic%20Field~TS05/ \nBuilt on PyGeopack https://github.com/mattkjames7/PyGeopack')
        def Bfield_n3(xvec, t):
            """Magnetic field model from Tsyganenko
            see self.model_name for specific model used
            """
            x, y, z = xvec.T
            time = to_utc(t)
            date_int, ut = time_to_tsyg(time)
            b = gp.ModelField(x, y, z, date_int, ut,
                                 Model=self.model_name,
                                 CoordIn=self.coord_in,
                                 CoordOut=self.coord_out)
            result = np.vstack(b).T
            return result

        @kamodofy(
                )
        def K_p(t):
            """K_p driver from Tsyganenko"""
            time = to_utc(t)
            date, ut = time_to_tsyg(time)
            params = gp.Params.GetModelParams(date, ut, Model = 'T89')
            return params['Kp'][0]
        
        @kamodofy(
                )
        def V(t):
            """Solar wind driver from Tsyganenko"""
            time = to_utc(t)
            date, ut = time_to_tsyg(time)
            params = gp.Params.GetModelParams(date, ut, Model = self.model_name)
            vx = params['Vx'][0]
            vy = params['Vy'][0]
            vz = params['Vz'][0]

            v = np.array([vx, vy, vz])
            return np.vstack(v).T
        
        
        in_expr = '\\vec{s}_{' + self.coord_in + '}'
        lhs_expr = '\\vec{r}_{' + self.coord_out + '} (' + in_expr + ')'
        soln_expr = lhs_expr + '=' + in_expr + ' + \\int_{s_0}^s \\hat{B} ( \\vec{r} (u)) du'
        
        @kamodofy(equation = '\\{ ' + lhs_expr + '_i \} \quad i \in [0,n-1], \quad ' + soln_expr,
                 hidden_args='n')
        def Bfield_trace(svec, t, n=100):

            time = to_utc(t)
            date_int, ut = time_to_tsyg(time)

            seeds = np.array(svec)
            if seeds.ndim == 1:  # If the array is 1D, it's assumed to have shape (3,)
                seeds = seeds.reshape(1, 3)

            for seed in seeds:
                x, y, z = seed
                T = self.trace(x, y, z, date_int, ut,
                               model_name=self.model_name,
                               coord_in=self.coord_in,
                               coord_out=self.coord_out)
                nans = np.isnan(T.s[0])
                
                x_pnts = get_component(T, 'x', self.coord_out.lower())[0][~nans]
                y_pnts = get_component(T, 'y', self.coord_out.lower())[0][~nans]
                z_pnts = get_component(T, 'z', self.coord_out.lower())[0][~nans]
                s = T.s[0][~nans]

                # make sure coords match self.coord_out
                pnts = np.vstack([x_pnts, y_pnts, z_pnts]).T
                
                s_default = np.linspace(s.min(), s.max(), n)
                
                @kamodofy(equation='\\begin{equation} ' + soln_expr + ' \\end{equation}')
                def f(t=s_default):
                    try:
                        return scipy.interpolate.interp1d(s, pnts, axis=0, fill_value='extrapolate')(t)
                    except ValueError as m:
                        print(m)
                        print(t, pnts[[0, -1]])
                        raise          
                    
                yield f

        
        
        self['Bvec'] = Bfield
        self['Bvec_n3'] = Bfield_n3
        self['F_B'] = Bfield_trace
        self['K_p'] = K_p
        self['V'] = V
# -

