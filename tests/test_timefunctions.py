import kamodo_ccmc.flythrough.model_wrapper as MW
import kamodo_ccmc.tools.timefunctions as tf
import datetime

ts = 1756771200.0
dt = datetime.datetime(2025, 9, 2, 0, 0, 0, tzinfo=datetime.timezone.utc)
iso = '2025-09-02T00:00:00.000000Z'
str = '2025/09/02 00:00:00'

def test00():
    '''This tests timeTStoSTR'''
    assert tf.timeTStoSTR(ts) == str

def test01():
    '''This tests timeTStoDT'''
    assert tf.timeTStoDT(ts) == dt

def test02():
    '''This tests timeDTtoISO'''
    assert tf.timeDTtoISO(dt) == iso

def test03():
    '''This tests timeISOtoDT'''
    assert tf.timeISOtoDT(iso) == dt

def test04():
    '''This tests timeDTtoSTR'''
    assert tf.timeDTtoSTR(dt) == str

def test05():
    '''This tests timeDTtoTS'''
    assert tf.timeDTtoTS(dt) == ts

def test06():
    '''This tests timeDTtoTS'''
    assert tf.timeDTtoTS(dt) == ts

