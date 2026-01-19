'''
Kamodo time conversion functions in this file:

 timeDTtoSTR
 timeDTtoISO
 timeDTtoTS
 timeISOtoDT
 timeTStoDT
 timeTStoSTR
 timeKOtoTS
 timeTStoKOoffset
'''

def timeDTtoSTR(inDT, datesep=None):
    '''
    UTC time conversion from datetime to formatted string
    
    IN:  2005-06-01 13:33:00.100000+00:00
    OUT: 2005/06/01 13:33:00
    '''
    
    STR = inDT.strftime("%Y/%m/%d %H:%M:%S")
    if datesep is not None:
        STR = STR.replace("/", datesep)
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
    
    IN:  2005-06-01T13:33:00.1Z  or  2005-152T13:33:00.1Z
    OUT: 2005-06-01 13:33:00.100000+00:00
    '''
    import dateutil.parser
    from datetime import datetime

    try:
        return dateutil.parser.parse(inISO)
    except (ValueError, dateutil.parser.ParserError):
        # Handle the YYYY-DDD format manually
        if '-' in inISO and 'T' in inISO:
            date_part, time_part = inISO.split('T')
            year, day_of_year = date_part.split('-')
            # Convert ordinal to standard date
            standard_date = datetime.strptime(f"{year}-{day_of_year}", "%Y-%j").strftime("%Y-%m-%d")
            return dateutil.parser.parse(f"{standard_date}T{time_part}")
        raise

def timeTStoDT(inTS):
    '''
    UTC time conversion from timestamp to datetime
    
    IN:  1117632780.1
    OUT: 2005-06-01 13:33:00.100000+00:00
    '''
    from datetime import datetime, timezone

    DT = datetime.fromtimestamp(float(inTS), tz=timezone.utc)

    return DT

def timeTStoSTR(inTS, datesep=None):
    '''
    UTC time conversion from timestamp to formatted string
    
    IN:  1117632780.1
    OUT: 2005/06/01 13:33:00
    '''
    from datetime import datetime, timezone

    DT = datetime.fromtimestamp(float(inTS), tz=timezone.utc)
    STR = DT.strftime("%Y/%m/%d %H:%M:%S")
    if datesep is not None:
        STR = STR.replace("/", datesep)
    return STR

def timeKOtoTS(ko, sOffset=0.):
    '''
    UTC time conversion from Kamodo object start (midnight) to timestamp
    
    IN:
      ko       A Kamodo object with ko.filedate defined (datetime)
      sOffset  Optional shift of timestamp by number of seconds given
    OUT:
      Timestamp, ie. 1117632780.1
    '''
    TS = ko.filedate.timestamp() + sOffset
    return TS

def timeTStoKOoffset(ko, inTS):
    '''
    UTC time conversion from timestamp to hours offset from start of Kamodo object

    IN:
      ko      A Kamodo object with ko.filedate defined (datetime)
      inTS    Timestamp in seconds to convert from
    OUT:
      returns floating point hours from start of Kamodo object
    '''
    Hrs = (float(inTS) - ko.filedate.timestamp()) / 3600.
    return Hrs

