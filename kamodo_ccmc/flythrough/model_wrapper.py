# -*- coding: utf-8 -*-
from datetime import datetime, timezone
from kamodo import get_defaults

model_dict = {'ADELPHI': 'AMPERE-Derived ELectrodynamic Properties of the ' +
                         'High-latitude Ionosphere ' +
                         'https://doi.org/10.1029/2020SW002677',
              'AMGeO': 'Assimilative Mapping of Geospace Observations ' +
                       'https://doi.org/10.5281/zenodo.3564914',
              'CTIPe': 'Coupled Thermosphere Ionosphere Plasmasphere ' +
                       'Electrodynamics Model ' +
                       'https://doi.org/10.1029/2007SW000364',
              'DTM': 'The Drag Temperature Model ' +
                     'https://doi.org/10.1051/swsc/2015001',
              'GAMERA_GM': 'Grid Agnostic MHD for Extended Research ' +
                        'Applications - Global Magnetosphere outputs ' +
                        'https://doi.org/10.3847/1538-4365/ab3a4c' +
                        ' (coming soon)',
              'GITM': 'Global Ionosphere Thermosphere Model ' +
                      'https://doi.org/10.1016/j.jastp.2006.01.008',
              'IRI': 'International Reference Ionosphere Model ' +
                     'https://doi.org/10.5194/ars-16-1-2018',
              'OpenGGCM_GM': 'The Open Geospace General Circulation Model - ' +
                             'Global Magnetosphere outputs only ' +
                             'https://doi.org/10.1023/A:1014228230714',
              'SuperDARN_uni': 'SuperDARN uniform grid output ' +
                               'https://doi.org/10.1029/2010JA016017',
              'SuperDARN_equ': 'SuperDARN equal area grid output ' +
                               'https://doi.org/10.1029/2010JA016017',
              'SWMF_IE': 'Space Weather Modeling Framework - Ionosphere and ' +
                         'Electrodynamics outputs ' +
                         'https://doi.org/10.1029/2006SW000272',
              'SWMF_GM': 'Space Weather Modeling Framework - Global ' +
                         'Magnetosphere outputs ' +
                         'https://doi.org/10.1029/2006SW000272',
              'TIEGCM': 'Thermosphere Ionosphere Electrodynamics General ' +
                        'Circulation Model ' +
                        'https://doi.org/10.1029/2012GM001297',
              'WACCMX': 'Whole Atmosphere Community Climate Model With ' +
                        'Thermosphere and Ionosphere Extension ' +
                        'https://doi.org/10.1002/2017MS001232',
              'WAMIPE': 'The coupled Whole Atmosphere Model - Ionosphere ' +
                        'Plasmasphere Model ' +
                        'https://doi.org/10.1002/2015GL067312 and ' +
                        'https://doi.org/10.1029/2022SW003193',
              'Weimer': 'Weimer Ionosphere model ' +
                        'https://doi.org/10.1029/2005JA011270'
              }


def Choose_Model(model=''):
    '''Returns module specific to the model requested.

    Input:
        model: A string corresponding to the desired model. If the
            string is empty, print statements are executed.
    Output: The module associated with the model reader for the desired model.
        If model is an empty string, None is returned.
    '''
    # UPDATE THIS AS MORE MODELS ARE ADDED

    if model == '':  # Give a list of possible values
        print("Possible models are: ")
        for key in model_dict.keys():
            print(f'{key}: {model_dict[key]}')
        return

    if model == 'CTIPe':
        import kamodo_ccmc.readers.ctipe_4D as module
        return module

    elif model == 'IRI':
        import kamodo_ccmc.readers.iri_4D as module
        return module

    elif model == 'GITM':
        import kamodo_ccmc.readers.gitm_4Dcdf as module
        return module

    elif model == 'SWMF_IE':
        import kamodo_ccmc.readers.swmfie_4Dcdf as module
        return module

    elif model == 'SWMF_GM':
        import kamodo_ccmc.readers.swmfgm_4D as module
        return module

    elif model == 'TIEGCM':
        import kamodo_ccmc.readers.tiegcm_4D as module
        return module

    elif model == 'OpenGGCM_GM':
        import kamodo_ccmc.readers.openggcm_gm_4Dcdf as module
        return module

    elif model == 'AMGeO':
        import kamodo_ccmc.readers.amgeo_4D as module
        return module

    elif model == 'SuperDARN_uni':
        import kamodo_ccmc.readers.superdarnuni_4D as module
        return module

    elif model == 'SuperDARN_equ':
        import kamodo_ccmc.readers.superdarnequ_4D as module
        return module

    elif model == 'ADELPHI':
        import kamodo_ccmc.readers.adelphi_4D as module
        return module

    elif model == 'WACCMX':
        import kamodo_ccmc.readers.waccmx_4D as module
        return module

    elif model == 'WAMIPE':
        import kamodo_ccmc.readers.wamipe_4D as module
        return module

    elif model == 'DTM':
        import kamodo_ccmc.readers.dtm_4D as module
        return module

    elif model == 'GAMERA_GM':
        import kamodo_ccmc.readers.gameragm_4D as module
        return module

    elif model == 'Weimer':
        import kamodo_ccmc.readers.weimer_4D as module
        return module

    else:
        raise AttributeError('Model not yet added: ' + str(model))


def Model_Reader(model):
    '''Returns model reader for requested model. Model agnostic.
    Input: model: A string or integer associated with the desired model.

    Output: The MODEL Kamodo class object in the desired model reader.
    '''

    module = Choose_Model(model)
    return module.MODEL()  # imports Kamodo


def Model_Variables(model, file_dir=None, return_dict=False):
    '''Returns model variables for requested model. Model agnostic.

    Inputs:
        model: A string or integer associated with the desired model.
        file_dir: A string giving the full file path for the chosen directory.
            The default value is None, meaning that tall of the variables
            possible for the chosen model will be returned/printed. Setting
            file_dir to a file path will restrict the out to only the variables
            found in the output files in the given directory.
        return_dict: A boolean (default=False). If False, the full variable
            dictionary associated with the model is printed to the screen and
            None is returned. If True, the full variable dictionary is returned
            without any printed statements.
        Output: None or the variable dictionary. (See return_dict description.)
    '''

    if file_dir == None:
        # choose the model-specific function to retrieve the variables
        module = Choose_Model(model)
        variable_dict = module.model_varnames
        var_dict = {value[0]: value[1:] for key, value in variable_dict.items()}
    
        # retrieve and print model specific and standardized variable names
        if return_dict:
            return var_dict
        else:
            print(f'\nThe {model} model accepts the standardized variable ' +
                  'names listed below.')
            print('-------------------------------------------------------------' +
                  '----------------------')
            for key, value in sorted(var_dict.items()):
                print(f"{key} : '{value}'")
            print()
            return
    else:
        reader = Model_Reader(model)
        ko = reader(file_dir, variables_requested='all')
    
        # either return or print nested_dictionary
        if return_dict:
            return ko.var_dict
        else:
            print('\nThe file directory contains the following ' +
                  'standardized variable names:')
            print('-------------------------------------------------' +
                  '----------------------------------')
            for key, value in ko.var_dict.items():
                print(f"{key} : '{value}'")
            return


def Variable_Search(search_string='', model='', file_dir='', return_dict=False):
    '''Search variable descriptions for the given string. If the model string
    is set, the chosen model will be searched. If file_dir is set to the
    directory where some model data is stored, the files available in that
    model data will be searched. If neither is set, then all model variables
    will be searched. Searching all of the model variables will take additional
    time as the model library grows. All search strings are converted to lower
    case text before search (e.g. "temperature" not "Temperature").
    
    Returns a dictionary with information that varies per call type if 
    return_dict is True. Default is to print the information to the screen
    and return nothing (return_dict=False).
    - If the search string, model and file directory are all not given or are
        all empty strings, then no dictionary is returned regardless of the
        value of the return_dict boolean.
    - if neither model nor file_dir is set, the search string is given, and
        the return_dict keyword is set to True, then the returned dictionary
        will have the model names as keys and the values will be dictionaries
        with the variable names as keys and the variable description,
        coordinate dependencies, and variable units as the value.
    - if only the model value is set, then the dictionary will have the
        variable names as keys and the variable description, coordinate
        dependencies, and variable units as the value.
    - if the model and file_dir values are set, then the dictionary will have
        the file names (or prefixes) as keys and the value will be a dictionary
        with the variable names as keys and the variable description,
        coordinate dependencies, and variable units as the value.
        '''

    search_string = search_string.lower()
    if search_string == '' and model == '' and file_dir == '':
        print('Printing all possible variables across all models...')
        for model in model_dict.keys():
            Model_Variables(model, return_dict=False)
        return None
    elif search_string == '' and model != '' and file_dir == '':
        new_dict = Model_Variables(model, return_dict=True)
        if not return_dict:
            for name, value in new_dict.items():
                print(name+':', value)
            return None
        else:
            return new_dict
    elif search_string == '' and model != '' and file_dir != '':
        new_dict = Model_Variables(model, file_dir=file_dir,
                                   return_dict=True)
        if not return_dict:
            for name, value in new_dict.items():
                print(name+':', value)
            return None
        else:
            return new_dict
    # remaining options assume search_string is set
    elif model == '' and file_dir == '':
        new_dict = {model: Variable_Search(search_string, model=model,
                                           return_dict=True) for 
                     model, desc in model_dict.items()}
        if not return_dict:
            for model, value in new_dict.items():
                print('\n'+model+':')
                if value == {}:
                    print(f'No {search_string} variables found.')
                else:
                    for name, values in value.items():
                        print(name+':', values)
            return None
        else:
            return new_dict
    elif file_dir == '' and model != '':
        var_dict = Model_Variables(model, return_dict=True)
        new_dict = {key: [value[0], value[-4]+'-'+value[-3], value[-2],
                          value[-1]] for key, value in
                    var_dict.items() if search_string in value[0].lower()}
        if not return_dict:
            for key, value in new_dict.items():
                print(key+':', value)
            return None
        else:
            return new_dict
    elif file_dir != '' and model != '':
        ko_var_dict = Model_Variables(model, file_dir, return_dict=True)
        new_dict = {name: [value[0], value[-4]+'-'+value[-3],
                           value[-2], value[-1]] for name, value in
                    ko_var_dict.items() if search_string in value[0].lower()}
        if new_dict == {}:
            print(f'No {search_string} variables found for {model} in ' +
                  f'{file_dir}.')
        if not return_dict:
            for name, value in new_dict.items():
                print(name+':', value)
            return None
        else:
            return new_dict
    else:
        return None


def File_Times(model, file_dir, print_output=True):
    '''Return/print the time range available in the data in the given dir. Also
    performs file conversions if the expected converted files are not found for
    the files or file naming patterns detected. If the model_times.txt and
    model_list.txt files are not found in the given directory, they are
    created. The user must delete these files to include any new files in the
    time range. Auto-detection of new files is no longer included to save
    execution cost in s3 buckets online.

    Inputs:
        model: A string or integer associated with the desired model.
        file_dir: A string giving the full file path for the chosen directory.
        print_output: A boolean (default=True). If True, the time range of the
            data in the given file directory is printed to the screen.
            If False, no print statements are executed.
    Output: Two datetime objects. The first datetime object is the date and
        time of the start time (UTC) of the beginning of the time grid, and
        the second is the same for the end of the time grid in the given
        file directory.
    '''

    # get time ranges from data
    from kamodo_ccmc.flythrough.SF_utilities import File_UTCTimes
    start_utcts, end_utcts, filedate = File_UTCTimes(model, file_dir)
    start_dt = datetime.utcfromtimestamp(start_utcts).replace(
        tzinfo=timezone.utc)
    end_dt = datetime.utcfromtimestamp(end_utcts).replace(tzinfo=timezone.utc)

    # print time ranges for given file groups: file_pattern, beg time, end time
    if print_output:
        print('UTC time ranges')
        print('------------------------------------------')
        print(f'Start {start_dt.strftime("Date: %Y-%m-%d  Time: %H:%M:%S")}')
        print(f'End {end_dt.strftime("Date: %Y-%m-%d  Time: %H:%M:%S")}')
    return start_dt, end_dt


def File_List(model, file_dir, print_output=False):
    '''Retrieve a list of the model output files in a given directory.
    Inputs:
        model: A string or integer associated with the desired model.
        file_dir: A string giving the full file path for the chosen directory.
        print_output: A boolean (default=True). If True, the time range of the
            data in the given file directory is printed to the screen.
            If False, no print statements are executed.
    Output: A list of the data files.
    '''

    reader = Model_Reader(model)
    kamodo_object = reader(file_dir, filetime=True, printfiles=print_output)
    return kamodo_object.filename


def Var_3D(model):
    '''Return list of model variables that are three-dimensional. Model
    agnostic.

    Input: model: A string or integer indicating the desired model.
    Output: A list of the standardized variable names associated with the model
        that have three dimensions (typically time + 2D spatial).
    '''

    # choose the model-specific function to retrieve the 3D variable list
    variable_dict = Model_Variables(model, return_dict=True)
    return [value[0] for key, value in variable_dict.items() if
            len(value[4]) == 3]


def Var_units(model, variable_list):
    '''Determines the proper units for the given variables.

    Inputs:
        model: A string indicating the desired model.
        variable_list: A list of strings for the desired standardized variable
            names.
    Output: A dictionary of the standardized variable names as keys with the
        corresponding units as the values.
    '''

    variable_dict = Model_Variables(model, return_dict=True)
    return {key: value[-1] for key, value in variable_dict.items() if key in
            variable_list}


def coord_units(coord_type, coord_grid):
    '''Determines the proper coordinate units given the coordinate system and
    type.

    Inputs:
        coord_type: A string corresponding to the desired
            coordinate system.
        coord_grid: A string corresponding to the desired
            coordinate grid type (e.g. 0 for cartesian or 1 for spherical.)
    Outputs: A dictionary with key, value pairs giving the proper units for
        each coordinate grid component (e.g. c1, c2, c3 in R_E for cartesian).
    '''

    if coord_grid == 'car':
        if coord_type in ['GDZ', 'SPH', 'RLL']:
            print(f'There is no cartesian version in the {coord_type} ' +
                  'coordinate system.')
            return
        return {'utc_time': 's', 'net_idx': '', 'c1': 'R_E', 'c2': 'R_E',
                'c3': 'R_E'}
    elif coord_grid == 'sph':
        if coord_type == 'GDZ':
            return {'utc_time': 's', 'net_idx': '', 'c1': 'deg', 'c2': 'deg',
                    'c3': 'km'}
        else:
            return {'utc_time': 's', 'net_idx': '', 'c1': 'deg', 'c2': 'deg',
                    'c3': 'R_E'}


def coord_names(coord_type, coord_grid):
    '''Determines the proper coordinate component names for a given coordinate
    system and grid type.

    Inputs:
        coord_type: An integer or string corresponding to the desired
            coordinate system.
        coord_grid: An integer or string corresponding to the desired
            coordinate grid type (e.g. 0 for cartesian or 1 for spherical.)
    Outputs: A dictionary with key, value pairs giving the proper coordinate
        component names for each coordinate grid component (e.g. c1, c2, c3
        associated with Longitude, Latitude, and Radius for most spherical
        systems).
    '''

    if coord_grid == 'car':
        return {'c1': 'X_'+coord_type, 'c2': 'Y_'+coord_type,
                'c3': 'Z_'+coord_type}
    elif coord_grid == 'sph' and coord_type == 'GDZ':
        return {'c1': 'Longitude', 'c2': 'Latitude', 'c3': 'Altitude'}
    elif coord_grid == 'sph' and coord_type != 'GDZ':
        return {'c1': 'Longitude', 'c2': 'Latitude', 'c3': 'Radius'}


def Coord_Range(kamodo_object, var_names, print_output=True,
                return_dict=False):
    '''Returns max and min of each dependent coordinate for the given variable
    in the given kamodo_object.

    Inputs:
        kamodo_object: an object returned by a kamodo model reader
        var_name: a list of variable names in the kamodo_object. For a complete
            list of variable names, type kamodo_object.detail().
        print_output: a boolean controlling whether information is printed to
            the screen, default is True.
        return_dict: a boolean controlling whether the constructed dictionary
            is returned, default is False.
    Output: A dictionary containing the max and min of each coordinate, named
        by the coordinate name.
    '''

    return_d = {}
    if print_output:
        print('The minimum and maximum values for each variable and ' +
              'coordinate are:', end="")
    for var in var_names:
        if print_output:
            print('\n' + var + ':')
        defaults = get_defaults(kamodo_object[var])
        key_list = list(defaults.keys())
        vals = defaults[key_list[0]].T
        if len(key_list) == 1 and len(vals.shape) > 1:  # not gridded
            if 'arg_units' in kamodo_object[var].meta.keys():  # coord names
                cunits = kamodo_object[var].meta['arg_units']
                names = list(cunits.keys())
                units = [value for key, value in cunits.items()]
                return_d[var] = {names[i]: [vals[i].min(), vals[i].max(),
                                            units[i]]
                                 for i in range(vals.shape[0])}
            else:
                return_d[var] = {'coord'+str(i): [vals[i].min(), vals[i].max()]
                                 for i in range(vals.shape[0])}
        else:  # assumes a gridded interpolator
            if 'arg_units' in kamodo_object[var].meta.keys():  # add units
                cunits = kamodo_object[var].meta['arg_units']
                return_d[var] = {key: [defaults[key].min(),
                                       defaults[key].max(),
                                       cunits[key.split('_')[0]]]
                                 for key in key_list}
                
            else:
                return_d[var] = {key: [defaults[key].min(),
                                       defaults[key].max()]
                                 for key in key_list}
        if print_output:
            for key in return_d[var]:
                print(key+':', return_d[var][key])
    if return_dict:
        return return_d
    else:
        return None
