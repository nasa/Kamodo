import kamodo_ccmc.flythrough.model_wrapper as MW
import kamodo_ccmc.flythrough.SatelliteFlythrough as SF
import types
import datetime
import pytest
from pathlib import Path
from math import isnan

model = 'CTIPe'
file_dir = 'TestData/'+model+'/'
variables_requested = ['T_n', 'T_e', 'TEC', 'E_theta300km']

def test00():
    '''
    This tests a list file can be found in output directory
    '''
    p = Path(file_dir+model+"_list.txt")
    assert p.is_file()

def test01_CTIPe_exists():
    '''
    This tests whether CTIPe exists as a model in kamodo
    '''
    assert type(MW.Choose_Model(model=model)) == types.ModuleType

def test02_CTIPe_variable():
    '''
    This tests whether a CTIPe variable search that includes "Temperature"
    has a variable "T_n" with units "K"
    '''
    vs = MW.Variable_Search('Temperature', model, return_dict=True)
    assert vs['T_n'][3] == 'K'

def test03_CTIPe_var_in_files():
    '''
    This tests that the variable "T_n" is in the test files
    '''
    vs = MW.Variable_Search('Temperature', model, file_dir, return_dict=True)
    assert vs['T_n'][3] == 'K'

def test04_CTIPe_times():
    '''
    This tests that proper start and end times are returned
    '''
    dt1 = datetime.datetime(2016, 1, 15, 0, 15, tzinfo=datetime.timezone.utc)
    dt2 = datetime.datetime(2016, 1, 17, 0, 0, tzinfo=datetime.timezone.utc)
    ft = MW.File_Times(model, file_dir)
    assert ft[0] == dt1 and ft[1] == dt2

def test05_CTIPe_interpolation():
    '''
    This tests creating a kamodo object, ko, and interpolating two different ways
    '''
    reader = MW.Model_Reader(model)
    ko = reader(file_dir, variables_requested=variables_requested[:1])
    if isnan(ko.T_n([5.2, 10., 60., 250.])[0]):
        raise AttributeError('Returned value is a NaN.')
    if isnan(ko.T_n_ijk(time=5.2, lon=10., lat=60., height=250.)):
        raise AttributeError('Returned value is a NaN.')
    if not ko.T_n([5.2, 10., 60., 250.]) == ko.T_n_ijk(time=5.2, lon=10., lat=60., height=250.):
        raise AttributeError('Values are not equal.')

def test06_CTIPe_coord_range():
    '''
    This tests coordinate range logic
    '''
    reader = MW.Model_Reader(model)
    ko = reader(file_dir)
    var_list = list(MW.Variable_Search('', model, file_dir, return_dict=True).keys())
    varijk_list = sorted(var_list + [item+'_ijk' for item in var_list])
    cr = MW.Coord_Range(ko, varijk_list, return_dict=True)
    assert cr['T_n']['time'][1] == pytest.approx(48.0, abs=.000001)

def test07_CTIPe_plot_value():
    '''
    This test makes a plotly figure and pulls a value out to compare to reference
    '''
    reader = MW.Model_Reader(model)
    ko = reader(file_dir, variables_requested=variables_requested)
    fig = ko.plot('T_n_ijk', plot_partial={'T_n_ijk': {'time': 40.0, 'height': 200.}})
    assert fig.data[0]['x'][2] == pytest.approx(-144.0, abs=.000001) and \
           fig.data[0]['y'][3] == pytest.approx(-84.0, abs=.000001) and \
           fig.data[0]['z'][4,5] == pytest.approx(822.16780419, abs=.000001)

def test08_CTIPe_flythrough():
    '''
    This tests simple flythrough extraction for a couple of points
    '''

