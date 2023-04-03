# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 10:50:21 2023

@author: rringuet
"""
# Adapted from:
# https://github.com/hapi-server/server-nodejs/blob/master/bin/Example.py

# Usage:
#   python2 Example.py --start START --stop STOP [-params PARAMETERS -fmt
#                                                 FORMAT]
#
# Generates a HAPI CSV file with a two parameters: a one-column scalar
# and a 3-column vector. The scalar values are the number of minutes
# since 1970-01-01. Vector component i at a given time is the scalar
# value at that time plus i. Only START and STOP are required and they
# must be in HAPI ISO8601 format.
#
# PARAMETERS is a comma-separated list of output parameters and
# can be any combination of Time, scalar, and vector. The default
# is 'Time,scalar,vector'. The order of parameters in the list is
# ignored and Time is always output as the first column.
#
# FORMAT is either 'csv' or 'binary'. The default is 'csv'.
#
# Examples:
#   python2 Example.py --start 1970-01-01T00:00:00.000000Z --stop
#                       1970-01-01T00:10:00.000000Z
#   python2 Example.py --fmt binary --start 1970-01-01T00:00:00.000000Z --stop
#                       1970-01-01T00:10:00.000000Z
#   python2 Example.py --params Time --start 1970-01-01T00:00:00.000000Z --stop
#                       1970-01-01T00:10:00.000000Z
#   python2 Example.py --params Time,vector --start 1970-01-01T00:00:00.000000Z
#                       --stop 1970-01-01T00:10:00.000000Z

import sys
import struct
import argparse
from datetime import datetime, timezone
import re
import urllib.request
import json
from numpy import vectorize, array, float32
import kamodo_ccmc.flythrough.model_wrapper as MW
from kamodo_ccmc.flythrough.SatelliteFlythrough import RealFlight

# On windows, error is
#   close failed in file object destructor:
#   sys.excepthook is missing
#   lost sys.stderr
# This is due to
#   https://bugs.python.org/issue11380
# The error message can be ignored.

# Trap broke pipe signal so usage in the form of
# python ./bin/Example.py | python lib/subset.py ...
# does not throw error when subset.py terminates read
# of output of Example.py.
if sys.platform != 'win32':
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

# get names from command given
parser = argparse.ArgumentParser()
parser.add_argument('--id', default='')
parser.add_argument('--params', default='')
parser.add_argument('--start', default='')
parser.add_argument('--stop', default='')
parser.add_argument('--fmt', default='csv')

# pull argument values. Dealing with times later.
v = vars(parser.parse_args())
# epoch = datetime.datetime(1970, 1, 1)  # not used
dataset = v['id']
params = v['params']
fmt = v["fmt"]

# get model, run name and satellite name from dataset
model_list = list(MW.model_dict.keys())  # list of model names
model_length = max([len(item) for item in model_list])  # max length of strings
with urllib.request.urlopen(
        "https://hapi-server.org/servers/SSCWeb/hapi/catalog") as url:
    data = json.load(url)
sat_names = [item['id'] for item in data['catalog']]  # list of sat names
# look for a model name at the beginning of the dataset name, with max length
model_mask = [model in dataset[:model_length] for model in model_list]
model = model_list[model_mask.index(True)]  # find first instance of model name
sat_name = dataset.split('_')[-1]   # sat name always the last piece
# run_name = dataset.split(model)[1].split(sat_name)[0][1:-1]
# model name will likely be in file_dir, so try below instead
file_dir = dataset[len(model)+1:].split(sat_name)[0][:-1]

# fake values for now, comment out before running for real testing
model, file_dir = 'CTIPe', 'D:/CTIPe/Data/'

# determine proper parameters list
# retrieve variable list and coord info for dataset
var_dict = MW.Variable_Search('', model=model, file_dir=file_dir,
                              return_dict=True)
coord_list = [value[4][-1] for key, value in var_dict.items()]
var_list = list(var_dict.keys())
# remove ilev variables from list b/c can't use in flythrough
if params == '':  # wants all variables
    variables = [var for i, var in enumerate(var_list) if coord_list[i]
                 not in ['ilev', 'ilev1']]
else:
    tmp = params.split(',')
    variables = [var for i, var in enumerate(var_list) if coord_list[i]
                 not in ['ilev', 'ilev1'] and var in tmp]
    if len(tmp) != len(variables):
        err_list = [var for var in tmp if var not in variables]
        print(f'The {err_list} variables were removed because they depend on' +
              'model-specific coordinate systems and cannot be used in ' +
              'the coordinate conversions required in the flythrough. Use ' +
              'the same variables names without the _ilev portion to get the' +
              ' desired variable.')
        if len(variables) == 0:
            raise AttributeError('No variables remain. Please choose a ' +
                                 'different list and try again.')
# remove spatial variables if requested and check for coordinate system
spatial_vars = [item for item in variables if 'X_' in item or 'Y_' in item
                or 'Z_' in item]
params_list = [item for item in variables if item not in spatial_vars]
if spatial_vars != []:
    coord_type = spatial_vars[0][2:]
else:
    coord_type = 'GEO'

# Start and stop times. Defaults are start and stop times of model output.
# RealFlight requires UTC timestamps as inputs.
start_dt, end_dt = MW.File_Times(model, file_dir, print_output=False)  # data
if v['start'] == '':
    start_utcts = start_dt.timestamp()
else:  # slice off 'Z' at the end
    start_utcts = datetime.fromisoformat(v['start'][:-1]).replace(
        tzinfo=timezone.utc).timestamp()
if v['stop'] == '':
    end_utcts = end_dt.timestamp()
else:
    end_utcts = datetime.fromisoformat(v['stop'][:-1]).replace(
        tzinfo=timezone.utc).timestamp()

# Run RealFlight to get data
results = RealFlight(sat_name, start_utcts, end_utcts, model, file_dir,
                     params_list, coord_type=coord_type)

# correct output variable list to include spatial variables if requested
if spatial_vars != []:
    strip_spatial = [item.split('_')[0] for item in spatial_vars]
    # replace spatial variables names with c1, c2, c3 as appropriate
    if 'X' in strip_spatial:
        strip_spatial[strip_spatial.index('X')] = 'c1'
    if 'Y' in strip_spatial:
        strip_spatial[strip_spatial.index('Y')] = 'c2'
    if 'Z' in strip_spatial:
        strip_spatial[strip_spatial.index('Z')] = 'c3'
    variables = strip_spatial + params_list
else:
    variables = params_list.copy()

# for now, all outputs are scalar because we separate out the components
'''
scalar = False
vector = False
if 'scalar' in params_list:
    scalar = True
if 'vector' in params_list:
    vector = True
'''

# convert UTC timestamps from results into iso strings for output
@vectorize
def ts_iso(utc_timestamp):
    '''Converts UTC timestamp from flythrough into datetime object.'''
    return datetime.utcfromtimestamp(utc_timestamp).isoformat()+'Z'


# rearrange output from RealFlight into an array for easier output
iso_timesZ = ts_iso(results['utc_time'])  # convert timestamps to iso strings
data = array([value for key, value in results.items() if key in variables],
             dtype=float32).T
del results  # clear memory because the results dict can be large

# set up output object and formatting string
if fmt == 'binary' and sys.version_info[0] == 3:  # Python 3
    import os
    stdout = os.fdopen(sys.stdout.fileno(), "wb", closefd=False)
    format_string = '<' + ''.join(['d' for item in data[0]])
elif fmt == 'csv':
    stdout = sys.stdout
    format_string = ''.join([',%d' for item in data[0]])
else:  # need to add JSON
    print('Formatting option not supported')

# perform output to given file
# NEED TO REPLACE NaNs WITH A FILL VALUE BEFORE OUTPUT!!!
for i in range(len(iso_timesZ)):
    # Python 3
    if fmt == 'binary':
        stdout.write(bytes("%s" % iso_timesZ[i], 'ascii'))
        stdout.write(struct.pack(format_string, *data[i]))
    elif fmt == 'csv':
        stdout.write("%s" % iso_timesZ[i])
        sys.stdout.write(format_string % tuple(data[i]))
        sys.stdout.write("\n")

if fmt == 'binary' and sys.version_info[0] == 3:
    stdout.flush()
