# HAPI_KamodoJSONs.py
"""
Created on Fri Mar 31 13:35:57 2023

@author: rringuet
"""
from numpy import double
import json
from kamodo import get_defaults
import kamodo_ccmc.flythrough.model_wrapper as MW

# user defined variables for JSON creation
dataset = 'CTIPe_D:/CTIPe/Storm_201303/'
variables_requested = ['v_nup_ilev', 'T_n', 'T_e', 'TEC', 'E_theta300km']
# Taken from the CTIPe model reader validation notebook with one extra.

# Demonstrate extraction of model and file_dir from dataset name.
# get list of model names in Kamodo
model_list = list(MW.model_dict.keys())
model_length = max([len(item) for item in model_list])
# retrieve model name from dataset
model_mask = [model in dataset[:model_length] for model in model_list]
model = model_list[model_mask.index(True)]
# retrieve run name from dataset
file_dir = dataset[len(model)+1:]
print(model, file_dir)

# Collect metadata from Kamodo
var_dict = MW.Variable_Search('', model, file_dir, return_dict=True)
var_list = list(var_dict.keys())
start_dt, stop_dt = MW.File_Times(model, file_dir, print_output=False)

# Create par_dict portion of json
# Likely need to add a fill value to replace the NaNs.
par_dict = [{'name': 'time', 'type': 'isotime', 'length': 24, 'units': 'UTC',
             'fill': 'null'}] +\
    [{'name': var, 'type': 'double', 'units': var_dict[var][-1],
      'description': var_dict[var][0] + ' in the '+var_dict[var][2] + '-' +
      var_dict[var][3] + ' coordinate system, dependent on ' +
      ''.join([item+', ' for item in var_dict[var][-2]])[:-2]}
     for var in var_list]

# Demonstrate creation of full dataset JSON and write to a file
json_dict = {
    "HAPI": "3.1",
    "status": {
        "code": 1200,
        "message": "OK"
    },
    "parameters": par_dict,
    "startDate": start_dt.isoformat().split('+')[0]+'Z',
    "stopDate": stop_dt.isoformat().split('+')[0]+'Z'}
with open(file_dir+model+'_dataset.json', 'w') as write_file:
    json.dump(json_dict, write_file)

# Demonstrate json creation for a given subset of variable names
new_par = [par_dict[0]] + [par_dict[i+1] for i, var in enumerate(var_list)
                           if var in variables_requested]
new_json = json_dict.copy()
new_json['parameters'] = new_par

# Still need to get the coordinate grids.
# Create the kamodo object for the dataset.
reader = MW.Model_Reader(model)
kamodo_object = reader(file_dir, variables_requested=variables_requested)

# Retrieve model spatial grids and store in the new json object
# reorder variables_requested to match json order for correct grid assignments
var_map = [var for var in var_list if var in variables_requested]
for i, var in enumerate(var_map):
    defaults = get_defaults(kamodo_object[var+'_ijk'])
    default_list = list(defaults.keys())
    new_json['parameters'][i+1]['size'] = [defaults[ll].shape[0] for ll in
                                           default_list[1:]]
    new_json['parameters'][i+1]['bins'] = \
        [{'name': key,
          'units': kamodo_object[var+'_ijk'].meta['arg_units'][key],
          'centers': list(double(defaults[key]))} for key in default_list[1:]
         ]

# Write parameter JSON
with open(file_dir+model+'_parameters.json', 'w') as write_file:
    json.dump(new_json, write_file)
print('JSON files created.')
