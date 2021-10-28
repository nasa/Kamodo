
import os

store = dict()


def initialize(fname, *variables):
    """Opens a kameleon-compatible file and initializes an interpolator"""

    # ccmc module must be imported here since it is unavailable in python3
    from ccmc import _CCMC as ccmc
    kameleon = ccmc.Kameleon()
    kameleon.open(fname)


    load_variables(kameleon, *variables)

    interpolator = kameleon.createNewInterpolator()

    store['interpolator'] = interpolator
    store['kameleon'] = kameleon

    metadata = get_global_metadata(kameleon)
    metadata['variables'] = get_variable_metadata(kameleon, *variables)
    return metadata

def get_global_metadata(kameleon):
    metadata = dict()
    for i in range(kameleon.getNumberOfGlobalAttributes()):
        gname = kameleon.getGlobalAttributeName(i)
        gattr = kameleon.getGlobalAttribute(gname)
        metadata[gname] = gattr.toString()
    return metadata

def get_variable_metadata(kameleon, *variables):
    metadata = dict()
    for varname in variables:
        metadata[varname] = dict()
        metadata[varname]['min'] = kameleon.getVariableAttribute(varname, 'actual_min').getAttributeFloat()
        metadata[varname]['max'] = kameleon.getVariableAttribute(varname, 'actual_max').getAttributeFloat()
        metadata[varname]['units'] = kameleon.getVisUnit(varname)
        metadata[varname]['native_units'] = kameleon.getNativeUnit(varname)
    return metadata

def create_interpolator():

    nvar = store['kameleon'].getNumberOfVariables()
    for i in range(nvar):
        varname = store['kameleon'].getVariableName(i)
        store['kameleon'].loadVariable(varname)
        
    store['interpolator'] = store['kameleon'].createNewInterpolator()
    return fname

def load_variables(kameleon, *variables):
    for var_name in variables:
        if kameleon.doesVariableExist(var_name):
            kameleon.loadVariable(var_name)         
        else:
            raise IOError('{} does not exist!'.format(var_name))

def interpolate(varname, *point):
    return store['interpolator'].interpolate(varname, *point)

if __name__ == '__channelexec__':
    for item in channel:
        try:
            channel.send(eval(item))
        except:
            if type(item) == tuple:
                channel.send(store['interpolator'].interpolate(
                    str(item[0]),
                    float(item[1]),
                    float(item[2]),
                    float(item[3])))
            elif type(item) == dict:
                if 'points' in item:
                    results = []
                    if 'variable' in item:
                        variable_name = str(item['variable'])
                        for point in item['points']:
                            results.append(store['interpolator'].interpolate(variable_name, *point))
                    elif 'variables' in item:
                        for point in item['points']:
                            result = []
                            for variable in item['variables']:
                                variable_name = str(variable)
                                result.append(store['interpolator'].interpolate(variable_name, *point))
                            results.append(result)
                    channel.send(results)
                else:
                    channel.send(item)
            else:
                channel.send(item)
