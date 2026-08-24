import kamodo_ccmc.flythrough.model_wrapper as MW
import types
import datetime
import pytest
from pathlib import Path
from math import isnan

# Anchor the path to the directory where this test script lives
TEST_DIR = Path(__file__).parent

model = 'SWMF_GM'
file_dir = str(TEST_DIR / 'TestData' / model) + '/'
variables_requested = ['B_z', 'B_y', 'B_x', 'theta_Btilt']

def test_data(get_test_data):
    '''
    This test makes sure data is downloaded
    '''
    model_dir = get_test_data(model)

def test_exists():
    '''
    This tests whether the model exists in kamodo
    '''
    assert type(MW.Choose_Model(model=model)) == types.ModuleType

def test_metadata():
    '''
    This tests metadata files can be recreated
    '''
    p1 = Path(file_dir+model+"_list.txt")
    p2 = Path(file_dir+model+"_times.txt")
    if p1.is_file(): p1.unlink()
    if p2.is_file(): p2.unlink()
    ft = MW.File_Times(model, file_dir)
    assert p1.is_file() and p2.is_file()

def test_times():
    '''
    This tests that proper start and end times are returned
    '''
    dt1 = datetime.datetime(2010, 12, 18, 0, 0, 0, tzinfo=datetime.timezone.utc)
    dt2 = datetime.datetime(2010, 12, 19, 23, 0, 0, tzinfo=datetime.timezone.utc)
    ft = MW.File_Times(model, file_dir)
    assert ft[0] == dt1 and ft[1] == dt2

def test_variable():
    '''
    This tests whether a variable search that includes "magnetic"
    has a variable "B_z" with units "nT"
    '''
    vs = MW.Variable_Search('magnetic', model, return_dict=True)
    assert vs['B_z'][3] == 'nT'

def test_var_in_files():
    '''
    This tests that the variable "B_z" is in the test files
    '''
    vs = MW.Variable_Search('magnetic', model, file_dir, return_dict=True)
    assert vs['B_z'][3] == 'nT'

def test_interpolation():
    '''
    This tests creating a kamodo object, ko, and interpolating two different ways
    '''
    reader = MW.Model_Reader(model)
    ko = reader(file_dir, variables_requested=variables_requested[:1])
    if isnan(ko.B_z([1.2, 10., 60., 50.])[0]):
        raise AttributeError('Returned value is a NaN.')
    if isnan(ko.B_z_ijk(time=1.2, X=10., Y=60., Z=50.)):
        raise AttributeError('Returned value is a NaN.')
    if not ko.B_z([1.2, 10., 60., 50.]) == ko.B_z_ijk(time=1.2, X=10., Y=60., Z=50.):
        raise AttributeError('Values are not equal.')

def test_coord_range():
    '''
    This tests coordinate range logic
    '''
    reader = MW.Model_Reader(model)
    ko = reader(file_dir)
    var_list = list(MW.Variable_Search('', model, file_dir, return_dict=True).keys())
    varijk_list = sorted(var_list + [item+'_ijk' for item in var_list])
    cr = MW.Coord_Range(ko, varijk_list, return_dict=True)
    assert cr['B_z']['time'][1] == 47.0

def test_plot_value():
    '''
    This test makes a plotly figure and pulls a value out to compare to reference
    '''
    reader = MW.Model_Reader(model)
    ko = reader(file_dir, variables_requested=variables_requested)
    fig = ko.plot('B_z_ijk', plot_partial={'B_z_ijk': {'time': 1.5, 'Z': 15.5}})
    assert fig.data[0]['x'][2] == pytest.approx(-204.0, abs=.000001) and \
           fig.data[0]['y'][3] == pytest.approx(-118.0, abs=.000001) and \
           fig.data[0]['z'][4,5] == pytest.approx(1.92658597, abs=.000001)

