import kamodo_ccmc.flythrough.model_wrapper as MW
import kamodo_ccmc.flythrough.SatelliteFlythrough as SF
import types
import datetime
import pytest
from pathlib import Path
from math import isnan

model = 'TIEGCM'
file_dir = 'TestData/'+model+'/'
variables_requested = ['T_n', 'rho', 'TEC', 'H_milev']

def test00():
    '''
    This tests a list file can be found in output directory
    '''
    p = Path(file_dir+model+"_list.txt")
    assert p.is_file()

def test01_exists():
    '''
    This tests whether the model exists in kamodo
    '''
    assert type(MW.Choose_Model(model=model)) == types.ModuleType

def test02_variable():
    '''
    This tests whether a variable search that includes "density"
    has a variable "rho" with units "g/cm**3"
    '''
    vs = MW.Variable_Search('density', model, return_dict=True)
    assert vs['rho'][3] == 'g/cm**3'

def test03_var_in_files():
    '''
    This tests that the variable "T_n" is in the test files
    '''
    vs = MW.Variable_Search('Temperature', model, file_dir, return_dict=True)
    assert vs['T_n'][3] == 'K'

def test04_times():
    '''
    This tests that proper start and end times are returned
    '''
    dt1 = datetime.datetime(2013, 3, 16, 0, 20, 0, 36, tzinfo=datetime.timezone.utc)
    dt2 = datetime.datetime(2013, 3, 17, 0, 0, tzinfo=datetime.timezone.utc)
    ft = MW.File_Times(model, file_dir)
    assert ft[0] == dt1 and ft[1] == dt2

def test05_interpolation():
    '''
    This tests creating a kamodo object, ko, and interpolating two different ways
    '''
    reader = MW.Model_Reader(model)
    ko = reader(file_dir, variables_requested=['T_n'])
    if isnan(ko.T_n([5.2, 10., 60., 350.])[0]):
        raise AttributeError('Returned value is a NaN.')
    if isnan(ko.T_n_ijk(time=5.2, lon=10., lat=60., height=350.)):
        raise AttributeError('Returned value is a NaN.')
    if not ko.T_n([5.2, 10., 60., 350.]) == ko.T_n_ijk(time=5.2, lon=10., lat=60., height=350.):
        raise AttributeError('Values are not equal.')

def test06_coord_range():
    '''
    This tests coordinate range logic
    '''
    reader = MW.Model_Reader(model)
    ko = reader(file_dir)
    var_list = list(MW.Variable_Search('', model, file_dir, return_dict=True).keys())
    varijk_list = sorted(var_list + [item+'_ijk' for item in var_list])
    cr = MW.Coord_Range(ko, varijk_list, return_dict=True)
    assert cr['T_n']['time'][1] == pytest.approx(24.0, abs=.000001)

def test07_plot_value():
    '''
    This test makes a plotly figure and pulls a value out to compare to reference
    '''
    reader = MW.Model_Reader(model)
    ko = reader(file_dir, variables_requested=['T_n'])
    fig = ko.plot('T_n_ijk', plot_partial={'T_n_ijk': {'time': 20.0, 'height': 200.}})
    assert fig.data[0]['x'][2] == pytest.approx(-170.0, abs=.000001) and \
           fig.data[0]['y'][3] == pytest.approx(-77.5, abs=.000001) and \
           fig.data[0]['z'][4,5] == pytest.approx(848.7894530288478, abs=.000001)

