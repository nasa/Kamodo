import inspect
from kamodo import Kamodo, kamodofy, gridify
from scipy.interpolate import RegularGridInterpolator

import cdflib
import pandas as pd
import numpy as np
import re

def time_interpolator(timekwarg):
    """{docstring}"""

    # Note: df will be passed into this function's local scope
    # t will be provisioned as a keyword argument

    df_ = df.reindex(df.index.union(t))
    df_interpolated = df_.interpolate(method='time')
    result = df_interpolated.reindex(t)
    return result

time_interpolator_docstring = """{varname} time interpolator

parameters:
t: datetime series or list

returns: {varname} [{units}] pandas DataFrame
"""

def get_interpolator(func, varname, data_frame, default_time, docstring):
    """Creates time interpolator with custom signature"""
    #extract source code from time_interpolator
    src = inspect.getsource(time_interpolator)

    #create variable-dependent signature
    new_src = (src \
           .format(docstring = docstring)
           .replace('timekwarg',"t = default_time")
           .replace('time_interpolator', varname))

#     default_time = self._instrument.data.index
    loc = dict(default_time = default_time, df = data_frame)
    exec(new_src, loc)
    return loc[varname]

def convert_ndimensional(data, index = None, columns = None):
    """converts high-dimensional data to a Dataframe"""
    if columns is None:
        columns = [range(i) for i in data.shape[1:]]
        columns = pd.MultiIndex.from_product(columns)

    return pd.DataFrame(data.T.reshape(data.shape[0], -1), 
        columns = columns, index = index)

class cdf_Kamodo(Kamodo):
    """Kamodofied cdflib

    Loading routines borrows heavily from pyspedas's cdf_to_tplot function
    """

    def __init__(self, filename, 
                varformat = '*', # regular expressions
                var_types = ['data', 'support_data'], 
                center_measurement = False,
                raise_errors = False,
                regnames = None,
                datetime = False,
                **kwargs):
        self._raise_errors = raise_errors
        self._filename = filename
        self._varformat = varformat
        self._var_types = var_types
        self._datetime = datetime
        self._var_types = var_types
        self._center_measurement = center_measurement

        #registration names map from file parameters to kamodo-compatible names
        if regnames is None:
            regnames = {}
        self._regnames = regnames 

        self._cdf_file = cdflib.CDF(self._filename)
        self._cdf_info = self._cdf_file.cdf_info()
        self.data = {} #python-in-Heliophysics Community data standard
        self.meta = {} #python-in-Heliophysics Community metadata standard
        self._dependencies = {}

        self._variable_names = self._cdf_info['rVariables'] +\
            self._cdf_info['zVariables']
        

        self._citation = self.get_citation()
        
        super(cdf_Kamodo, self).__init__(**kwargs)
        
        self.load_variables()
        self.register_variables()
        
        
    def get_dependency(self, x_axis_var):
        """Retrieves variable dependency unique to filename"""
        return self._dependencies.get(self._filename + x_axis_var)

    def set_dependency(self, x_axis_var, x_axis_data):
        """Sets variable dependency unique to filename"""
        self._dependencies[self._filename + x_axis_var] = x_axis_data
    
    def set_epoch(self, x_axis_var):
        """Stores epoch dependency"""

        data_type_description \
            = self._cdf_file.varinq(x_axis_var)['Data_Type_Description']

        center_measurement = self._center_measurement
        cdf_file = self._cdf_file
        if self.get_dependency(x_axis_var) is None:
            delta_plus_var = 0.0
            delta_minus_var = 0.0
            delta_time = 0.0

            xdata = cdf_file.varget(x_axis_var)
            epoch_var_atts = cdf_file.varattsget(x_axis_var)

            # check for DELTA_PLUS_VAR/DELTA_MINUS_VAR attributes
            if center_measurement:
                if 'DELTA_PLUS_VAR' in epoch_var_atts:
                    delta_plus_var = cdf_file.varget(
                        epoch_var_atts['DELTA_PLUS_VAR'])
                    delta_plus_var_att = cdf_file.varattsget(
                        epoch_var_atts['DELTA_PLUS_VAR'])

                    # check if a conversion to seconds is required
                    if 'SI_CONVERSION' in delta_plus_var_att:
                        si_conv = delta_plus_var_att['SI_CONVERSION']
                        delta_plus_var = delta_plus_var.astype(float) \
                            * np.float(si_conv.split('>')[0])
                    elif 'SI_CONV' in delta_plus_var_att:
                        si_conv = delta_plus_var_att['SI_CONV']
                        delta_plus_var = delta_plus_var.astype(float) \
                            * np.float(si_conv.split('>')[0])

                if 'DELTA_MINUS_VAR' in epoch_var_atts:
                    delta_minus_var = cdf_file.varget(
                        epoch_var_atts['DELTA_MINUS_VAR'])
                    delta_minus_var_att = cdf_file.varattsget(
                        epoch_var_atts['DELTA_MINUS_VAR'])

                    # check if a conversion to seconds is required
                    if 'SI_CONVERSION' in delta_minus_var_att:
                        si_conv = delta_minus_var_att['SI_CONVERSION']
                        delta_minus_var = \
                            delta_minus_var.astype(float) \
                            * np.float(si_conv.split('>')[0])
                    elif 'SI_CONV' in delta_minus_var_att:
                        si_conv = delta_minus_var_att['SI_CONV']
                        delta_minus_var = \
                            delta_minus_var.astype(float) \
                            * np.float(si_conv.split('>')[0])

                # sometimes these are specified as arrays
                if isinstance(delta_plus_var, np.ndarray) \
                        and isinstance(delta_minus_var, np.ndarray):
                    delta_time = (delta_plus_var
                                  - delta_minus_var) / 2.0
                else:  # and sometimes constants
                    if delta_plus_var != 0.0 or delta_minus_var != 0.0:
                        delta_time = (delta_plus_var
                                      - delta_minus_var) / 2.0

        if self.get_dependency(x_axis_var) is None:
            if ('CDF_TIME' in data_type_description) or \
                    ('CDF_EPOCH' in data_type_description):
                xdata = cdflib.cdfepoch.unixtime(xdata)
                xdata = np.array(xdata) + delta_time
                if self._datetime:
                    xdata = pd.to_datetime(xdata,  unit = 's')
                self.set_dependency(x_axis_var, xdata)

    def get_index(self, variable_name):
        var_atts = self._cdf_file.varattsget(variable_name)

        if "DEPEND_TIME" in var_atts:
            x_axis_var = var_atts["DEPEND_TIME"]
            self.set_epoch(x_axis_var)
        elif "DEPEND_0" in var_atts:
            x_axis_var = var_atts["DEPEND_0"]
            self.set_epoch(x_axis_var)

        dependencies = []
        for suffix in ['TIME'] + list('0123'):
            dependency = "DEPEND_{}".format(suffix)
            dependency_name = var_atts.get(dependency)
            if dependency_name is not None:
                dependency_data = self.get_dependency(dependency_name)
                if dependency_data is None:
                    dependency_data = self._cdf_file.varget(dependency_name)
                    # get first unique row
                    dependency_data = pd.DataFrame(dependency_data).drop_duplicates().values[0]
                    self.set_dependency(dependency_name, dependency_data)
                dependencies.append(dependency_data)

        index_ = None
        if len(dependencies) == 0:
            pass
        elif len(dependencies) == 1:
            index_ = dependencies[0]
        else:
            index_ = pd.MultiIndex.from_product(dependencies)

        return index_


    def load_variables(self):
        """loads cdf variables based on varformat

        Based heavily on cdf_to_tplot from pyspedas
        """
        varformat = self._varformat
        if varformat is None:
            varformat = ".*"

        varformat = varformat.replace("*", ".*")
        var_regex = re.compile(varformat)

        for variable_name in self._variable_names:
            if not re.match(var_regex, variable_name):
                # skip this variable
                continue
            var_atts = self._cdf_file.varattsget(variable_name)

            if 'VAR_TYPE' not in var_atts:
#                 print('skipping {} (no VAR_TYPE)'.format(variable_name))
                continue

            if var_atts['VAR_TYPE'] not in self._var_types:
#                 print('skipping {} ({})'.format(variable_name, var_atts['VAR_TYPE']))
                continue

            var_properties = self._cdf_file.varinq(variable_name)

            try:
                ydata = self._cdf_file.varget(variable_name)
            except (TypeError):
#                 print('skipping {} (TypeError)'.format(variable_name))
                continue

            if ydata is None:
#                 print('skipping {} (empty)'.format(variable_name))
                continue
                

            if "FILLVAL" in var_atts:
                if (var_properties['Data_Type_Description'] == 'CDF_FLOAT'
                    or var_properties['Data_Type_Description']
                    == 'CDF_REAL4'
                    or var_properties['Data_Type_Description']
                    == 'CDF_DOUBLE'
                    or var_properties['Data_Type_Description']
                    == 'CDF_REAL8'):
                    if ydata[ydata == var_atts["FILLVAL"]].size != 0:
                        ydata[ydata == var_atts["FILLVAL"]] = np.nan
            
            index = self.get_index(variable_name)

            
            try:
                if isinstance(index, pd.MultiIndex):
                    self.data[variable_name] = pd.DataFrame(ydata.ravel(), index = index)
                else:
                    if len(ydata.shape) == 1:
                        self.data[variable_name] = pd.Series(ydata, index = index)
                    elif len(ydata.shape) == 2:
                        self.data[variable_name] = pd.DataFrame(ydata, index = index)
                    elif len(ydata.shape) >2:
                        self.data[variable_name] = convert_ndimensional(ydata, index = index)
                    else:
                        raise NotImplementedError('Cannot handle {} with shape {}'.format(variable_name, ydata.shape))
            except:
                self.data[variable_name] = {'ydata':ydata, 'index':index}
                if self._raise_errors:
                    raise

            self.meta[variable_name] = var_atts
            

    def get_citation(self, docs = ['Project', 
                                    'Source_name', 
                                    'PI_name', 
                                    'PI_affiliation',  
                                    'TEXT',
                                    'Rules_of_use',
                                    'HTTP_LINK']):
        """Extracts citation info from metadata"""

        global_atts = self._cdf_file.globalattsget()
        citation = ''
        for k in docs:
            v = global_atts.get(k, '')
            if type(v) is str:
                vstr = v
            else:
                vstr = '\n\t'.join(v)

            citation += "{}:\n\t{}\n".format(k, vstr)
        return citation

    def register_variables(self):
        for variable_name, df in self.data.items():
        
            dependencies = []
            for i in ['TIME'] + list('01234'):
                dependency_name = self.meta[variable_name].get('DEPEND_{}'.format(i))
                if dependency_name is not None:
                    dependencies.append(dependency_name)
            if len(dependencies) == 0:
                # print('not registering {}: no dependencies'.format(variable_name))
                continue

            if not hasattr(df, 'index'):
                print('{} has no index, skipping..'.format(variable_name))
                continue

            regname = self._regnames.get(variable_name, variable_name)
            docstring = self.meta[variable_name]['CATDESC']
            units = self.meta[variable_name]['UNITS']
            citation = self._citation


            if isinstance(df.index, pd.MultiIndex):
                indices = df.index.levels
                data_shape = [len(i) for i in indices]
                grid_interpolator = RegularGridInterpolator(
                    indices, 
                    df.values.reshape(data_shape),
                    bounds_error = False)

                grid_interpolator.__name__ = variable_name

                grid_args = {d: self.get_dependency(d) for d in dependencies}
                interpolator = gridify(grid_interpolator, **grid_args)

            # catch time dependencies such as Epoch_state and Epoch
            elif (len(dependencies) == 1) & ('epoch' in dependencies[0].lower()):
                interpolator = get_interpolator(time_interpolator, 
                                                regname,
                                                df,
                                                df.index,
                                                docstring)
            else:
                print('can not register {}'.format(variable_name))
                continue


            
            self[regname] = kamodofy(interpolator, 
                                     units = units,
                                     citation = citation)
    
                