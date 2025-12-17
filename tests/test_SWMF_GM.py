import kamodo_ccmc.flythrough.model_wrapper as MW
import types
model = 'SWMF_GM'

def test00_SWMF_GM_exists():
    '''This tests whether SWMF_GM exists as a model in kamodo'''
    assert type(MW.Choose_Model(model=model)) == types.ModuleType

def test01_SWMF_GM_variable():
    '''This tests whether a SWMF_GM variable search that includes "magnetic" has a variable "B_z" with units "nT"'''
    assert MW.Variable_Search('magnetic', model, return_dict=True)['B_z'][3] == 'nT'