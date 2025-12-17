import kamodo_ccmc.flythrough.model_wrapper as MW
import types
model = 'GITM'

def test00_GITM_exists():
    '''This tests whether GITM exists as a model in kamodo'''
    assert type(MW.Choose_Model(model=model)) == types.ModuleType

def test01_GITM_variable():
    '''This tests whether a GITM variable search that includes "Temperature" has a variable "T_n" with units "K"'''
    assert MW.Variable_Search('Temperature', model, return_dict=True)['T_n'][3] == 'K'
