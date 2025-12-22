import kamodo_ccmc.flythrough.model_wrapper as MW
import kamodo_ccmc.flythrough.SatelliteFlythrough as SF
import types
import datetime
import pytest
from pathlib import Path
from math import isnan

model = 'AMGeO'
file_dir = 'TestData/'+model+'/'
variables_requested = ['Sigma_H', 'W_JouleN', 'W_JouleS', 'B_z']

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
    This tests whether a variable search that includes "Joule"
    has a variable "Q_Joule" with units "mW/m**2"
    '''
    vs = MW.Variable_Search('Joule', model, return_dict=True)
    assert vs['Q_Joule'][3] == 'mW/m**2'

def test03_var_in_files():
    '''
    This tests that the variable "Q_Joule" is in the test files
    '''
    vs = MW.Variable_Search('Joule', model, file_dir, return_dict=True)
    assert vs['Q_Joule'][3] == 'mW/m**2'

def test04_times():
    '''
    This tests that proper start and end times are returned
    '''
    dt1 = datetime.datetime(2015, 6, 22, 0, 2, 0, 6, tzinfo=datetime.timezone.utc)
    dt2 = datetime.datetime(2015, 6, 23, 23, 57, 0, 2747, tzinfo=datetime.timezone.utc)
    ft = MW.File_Times(model, file_dir)
    assert ft[0] == dt1 and ft[1] == dt2

def test05_interpolation():
    '''
    This tests creating a kamodo object, ko, and interpolating two different ways
    '''
    reader = MW.Model_Reader(model)
    ko = reader(file_dir, variables_requested=variables_requested[:1])
    if isnan(ko.Sigma_H([5.2, 10., 60.])[0]):
        raise AttributeError('Returned value is a NaN.')
    if isnan(ko.Sigma_H_ijk(time=5.2, lon=10., lat=60.)):
        raise AttributeError('Returned value is a NaN.')
    if not ko.Sigma_H([5.2, 10., 60.]) == ko.Sigma_H_ijk(time=5.2, lon=10., lat=60.):
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
    assert cr['Q_Joule']['time'][1] == pytest.approx(47.95000076293945, abs=.000001)

def test07_plot_value():
    '''
    This test makes a plotly figure and pulls a value out to compare to reference
    '''
    reader = MW.Model_Reader(model)
    ko = reader(file_dir, variables_requested=variables_requested)
    fig = ko.plot('Sigma_H_ijk', plot_partial={'Sigma_H_ijk': {'lat': 60.0}})
    assert fig.data[0]['x'][2] == pytest.approx(0.2, abs=.000001) and \
           fig.data[0]['y'][3] == pytest.approx(-150.0, abs=.000001) and \
           fig.data[0]['z'][4,5] == pytest.approx(5.745647529104042, abs=.000001)

