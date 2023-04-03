# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 13:35:57 2023

@author: rringuet
"""
import json
import urllib.request
import kamodo_ccmc.flythrough.model_wrapper as MW

# user determined inputs for JSON creation
dataset = 'CTIPe_D:/CTIPe/Storm_201303/_grace1'
variables_requested = ['v_nup_ilev', 'T_n', 'T_e', 'TEC', 'E_theta300km']
# Taken from the CTIPe model reader validation notebook with one extra f(ilev).

# Demonstrate extraction of model and file_dir from dataset name.
# get list of model names in Kamodo
model_list = list(MW.model_dict.keys())
model_length = max([len(item) for item in model_list])
# retrieve model name from dataset
model_mask = [model in dataset[:model_length] for model in model_list]
model = model_list[model_mask.index(True)]
# get list of satellite names from SSCWeb
with urllib.request.urlopen(
        "https://hapi-server.org/servers/SSCWeb/hapi/catalog") as url:
    data = json.load(url)
sat_names = [item['id'] for item in data['catalog']]
# split dataset into model, file_dir, and sat_name
sat_name = dataset.split('_')[-1]
# retrieve run name from dataset
file_dir = dataset[len(model)+1:].split(sat_name)[0][:-1]
print(model, file_dir, sat_name)

# Collect metadata from Kamodo
var_dict = MW.Variable_Search('', model, file_dir, return_dict=True)
var_list = list(var_dict.keys())
start_dt, stop_dt = MW.File_Times(model, file_dir, print_output=False)

# Create par_dict portion of json
# Likely need to add a fill value to replace the NaNs.
par_dict = [{'name': 'time', 'type': 'isotime', 'length': 24, 'units': 'UTC',
             'fill': 'null'}] +\
    [{'name': var, 'type': 'double', 'units': var_dict[var][-1],
      'description': var_dict[var][0] + ' in the ' + var_dict[var][2] + '-' +
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

# Demonstrate json creation for a given subset of variable names and write
new_par = [par_dict[0]] + [par_dict[i+1] for i, var in enumerate(var_list)
                           if var in variables_requested]
new_json = json_dict.copy()
new_json['parameters'] = new_par
with open(file_dir+model+'_parametersRR.json', 'w') as write_file:
    json.dump(new_json, write_file)
print('JSON files created.')
