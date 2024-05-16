def timeDTtoSTR(inDT):
    '''
    UTC time conversion from datetime to formatted string
    
    IN:  2005-06-01 13:33:00.100000+00:00
    OUT: 2005/06/01 13:33:00
    '''
    
    STR = inDT.strftime("%Y/%m/%d %H:%M:%S")
    return STR

def timeDTtoISO(inDT):
    '''
    UTC time conversion from datetime to ISO
    
    IN:  2005-06-01 13:33:00.100000+00:00
    OUT: 2005-06-01T13:33:00.1Z
    '''
    
    ISO = inDT.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    return ISO

def timeDTtoTS(inDT):
    '''
    UTC time conversion from datetime to timestamp
    
    IN:  2005-06-01 13:33:00.100000+00:00
    OUT: 1117632780.1
    '''

    TS = inDT.timestamp()
    return TS

def timeISOtoDT(inISO):
    '''
    UTC time conversion from ISO string to datetime
    
    IN:  2005-06-01T13:33:00.1Z
    OUT: 2005-06-01 13:33:00.100000+00:00
    '''
    from dateutil import parser

    DT = parser.parse(inISO)
    return DT

def timeTStoDT(inTS):
    '''
    UTC time conversion from timestamp to datetime
    
    IN:  1117632780.1
    OUT: 2005-06-01 13:33:00.100000+00:00
    '''
    import pytz
    from datetime import datetime

    DT = pytz.utc.localize(datetime.utcfromtimestamp(inTS))
    return DT

def timeTStoSTR(inTS):
    '''
    UTC time conversion from timestamp to formatted string
    
    IN:  1117632780.1
    OUT: 2005/06/01 13:33:00
    '''
    import pytz
    from datetime import datetime

    DT = pytz.utc.localize(datetime.utcfromtimestamp(inTS))
    STR = DT.strftime("%Y/%m/%d %H:%M:%S")
    return STR

def timeKOtoTS(ko, sOffset=0.):
    '''
    UTC time conversion from Kamodo object start (midnight) to timestamp
    
    IN:      A Kamodo object with ko.filedate defined (datetime)
    sOffset: Optional shift of timestamp by number of seconds given
    OUT:     1117632780.1
    '''
    TS = ko.filedate.timestamp() + sOffset
    return TS

