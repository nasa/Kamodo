#!/usr/bin/env python3

#
# NOSA HEADER START
#
# The contents of this file are subject to the terms of the NASA Open
# Source Agreement (NOSA), Version 1.3 only (the "Agreement").  You may
# not use this file except in compliance with the Agreement.
#
# You can obtain a copy of the agreement at
#   docs/NASA_Open_Source_Agreement_1.3.txt
# or
#   https://sscweb.gsfc.nasa.gov/WebServices/NASA_Open_Source_Agreement_1.3.txt.
#
# See the Agreement for the specific language governing permissions
# and limitations under the Agreement.
#
# When distributing Covered Code, include this NOSA HEADER in each
# file and include the Agreement file at
# docs/NASA_Open_Source_Agreement_1.3.txt.  If applicable, add the
# following below this NOSA HEADER, with the fields enclosed by
# brackets "[]" replaced with your own identifying information:
# Portions Copyright [yyyy] [name of copyright owner]
#
# NOSA HEADER END
#
# Copyright (c) 2013-2019 United States Government as represented by
# the National Aeronautics and Space Administration. No copyright is
# claimed in the United States under Title 17, U.S.Code. All Other
# Rights Reserved.
#

"""
Example of accessing the Satellite Situation Center (SSC) web services
https://sscweb.gsfc.nasa.gov/WebServices/REST/.
"""

import os
import math
import platform
import xml.etree.ElementTree as ET
from datetime import datetime, timezone, timedelta
from enum import Enum
#import logging
import httplib2
import dateutil.parser
import numpy as np


__version__="0.1.1"


class CoordinateSystem(Enum):
    """
    Python Enum representing the CoordinateSystem type defined
    in SSC.xsd.
    """
    GEO = 'Geo'
    GM = 'Gm'
    GSE = 'Gse'
    GSM = 'Gsm'
    SM = 'Sm'
    GEI_TOD = 'GeiTod'
    GEI_J_2000 = 'GeiJ2000'


class CoordinateSystemType(Enum):
    """
    Python Enum representing the CoordinateSystemType defined in
    SSC.xsd.
    """
    SPHERICAL = 'Spherical'
    CARTESIAN = 'Cartesian'


class CoordinateComponent(Enum):
    """
    Python Enum representing the CoordinateComponent type defined in
    SSC.xsd.
    """
    X = 'X'
    Y = 'Y'
    Z = 'Z'
    LAT = 'Lat'
    LON = 'Lon'
    LOCAL_TIME = 'Local_Time'


class InternalBFieldModel(Enum):
    """
    Python Enum representing the InteralBFieldModel type defined in
    SSC.xsd.
    """
    IGRF = 'IGRF'
    SIMPLE_DIPOLE = 'SimpleDipole'


class ExternalBFieldModel(Enum):
    """
    Python Enum representing the ExternalBFieldModel type defined in
    SSC.xsd.
    """
    TSYGANENKO96 = 'Tsyganenko96BFieldModel'
    TSYGANENKO89C = 'Tsyganenko89cBFieldModel'
    TSYGANENKO87 = 'Tsyganenko87BFieldModel'


class Tsyganenko87Kp(Enum):
    """
    Python Enum representing the Tsyganenko87Kp type defined in
    SSC.xsd.
    """
    KP_0_0 = 'KP0_0'
    KP_1_1_1 = 'KP1_1_1'
    KP_2_2_2 = 'KP2_2_2'
    KP_3_3_3 = 'KP3_3_3'
    KP_4_4_4 = 'KP4_4_4'
    KP_5 = 'KP5'


class Tsyganenko89cKp(Enum):
    """
    Python Enum representing the Tsyganenko89cKp type defined in
    SSC.xsd.
    """
    KP_0_0 = 'KP0_0'
    KP_1_1_1 = 'KP1_1_1'
    KP_2_2_2 = 'KP2_2_2'
    KP_3_3_3 = 'KP3_3_3'
    KP_4_4_4 = 'KP4_4_4'
    KP_5_5_5 = 'KP5_5_5'
    KP_6 = 'KP6'


class BFieldTraceDirection(Enum):
    """
    Python Enum representing the BFieldTraceDirection type defined
    in SSC.xsd.
    """
    SAME_HEMISPHERE = 'SameHemisphere'
    OPPOSITE_HEMISPHERE = 'OppositeHemisphere'
    NORTH_HEMISPHERE = 'NorthHemisphere'
    SOUTH_HEMISPHERE = 'SouthHemisphere'
    EITHER_HEMISPHERE = 'EitherHemisphere'


class ConditionOperator(Enum):
    """
    Python Enum representing the ConditionOperator type defined in
    SSC.xsd.
    """
    ALL = 'All'
    ANY = 'Any'


class ConjunctionAreaType(Enum):
    """
    Python Enum representing the ConjunctionAreaType defined in SSC.xsd.
    """
    GEO_BOX = 'GeoBox'
    GM_BOX = 'GmBox'
    DISTANCE = 'Distance'


class DateFormat(Enum):
    """
    Python Enum representing the DateFormat type defined in SSC.xsd.
    """
    YYYY_DDD = 'yyyy_ddd'
    YY_MM_DD = 'yy_mm_dd'
    YY_MMM_DD = 'yy_Mmm_dd'
    YY_CMMM_DD = 'yy_CMMM_dd'


class DegreeFormat(Enum):
    """
    Python Enum representing the DegreeFormat type defined in SSC.xsd.
    """
    DECIMAL = 'Decimal'
    MINUTES = 'Minutes'
    MINUTES_SECONDS = 'MinutesSeconds'

class DistanceFormat(Enum):
    """
    Python Enum representing the DistanceFormat type defined in SSC.xsd.
    """
    RE = 'Re'
    KM = 'Km'
    INTEGER_KM = 'IntegerKm'
    SCIENTFIC_NOTATION_KM = 'ScientificNotationKm'

class FootpointRegion(Enum):
    """
    Python Enum representing the FootpointRegion type defined in SSC.xsd.
    """
    NOT_APPLICABLE = 'NotApplicable'
    NORTH_CUSP = 'NorthCusp'
    SOUTH_CUSP = 'SouthCusp'
    NORTH_CLEFT = 'NorthCleft'
    SOUTH_CLEFT = 'SouthCleft'
    NORTH_AURORAL_OVAL = 'NorthAuroralOval'
    SOUTH_AURORAL_OVAL = 'SouthAuroralOval'
    NORTH_POLAR_CAP = 'NorthPolarCap'
    SOUTH_POLAR_CAP = 'SouthPolarCap'
    NORTH_MID_LATITUDE = 'NorthMidLatitude'
    SOUTH_MID_LATITUDE = 'SouthMidLatitude'
    LOW_LATITUDE = 'LowLatitude'

class Hemisphere(Enum):
    """
    Python Enum representing the Hemisphere type defined in SSC.xsd.
    """
    SOUTH = 'South'
    NORTH = 'North'

class LatLonFormat(Enum):
    """
    Python Enum representing the LatLonFormat type defined in SSC.xsd.
    """
    LAT_90_LON_360 = 'Lat90Lon360'
    LAT_90_LON_180 = 'Lat90Lon180'
    LAT_90_SN_LON_180_WE = 'Lat90SnLon180We'

class MapProjection(Enum):
    """
    Python Enum representing the MapProjection type defined in SSC.xsd.
    """
    AZIMUTHAL = 'Azimuthal'
    CYLINDRICAL = 'Cylindrical'
    MERCATOR = 'Mercator'
    MOLLEWEIDE = 'Molleweide'
    ORTHOGRAPHIC = 'Othographic'
    STEREOGRAPHIC = 'Stereographic'

class MapProjectionTrace(Enum):
    """
    Python Enum representing the MapProjectionTrace type defined
    in SSC.xsd.
    """
    B_FIELD_NORTH = 'BFieldNorth'
    B_FIELD_SOUTH = 'BFieldSouth'
    RADIAL = 'Radial'

class MapRegion(Enum):
    """
    Python Enum representing the MapRegion type defined in SSC.xsd.
    """
    NORTH_CUSP = 'NorthCusp'
    SOUTH_CUSP = 'SouthCusp'
    NORTH_CLEFT = 'NorthCleft'
    SOUTH_CLEFT = 'SouthCleft'
    NORTH_AURORAL_OVAL = 'NorthAuroralOval'
    SOUTH_AURORAL_OVAL = 'SouthOuroralOval'
    NORTH_POLAR_CAP = 'NorthPolarCap'
    SOUTH_POLAR_CAP = 'SouthPolarCap'
    NORTH_MID_LATITUDE = 'NorthMidLatitude'
    SOUTH_MID_LATITUDE = 'SouthMidLatitude'
    LOW_LATITUDE = 'LowLatitude'
    NONE = 'None'

class PolarMapOrientation(Enum):
    """
    Python Enum representing the PolarMapOrientation type defined in
    SSC.xsd.
    """
    EQUATORIAL = 'Equatorial'
    NORTH_POLE = 'NorthPole'
    SOUTH_POLE = 'SouthPole'

class ProjectionCoordinateSystem(Enum):
    """
    Python Enum representing the ProjectionCoordinateSystem type defined
    in SSC.xsd.
    """
    GEO = 'Geo'
    GM = 'Gm'
    SM = 'Sm'

class QueryResultType(Enum):
    """
    Python Enum representing the QueryResultType defined in SSC.xsd.
    """
    XML = 'Xml'
    LISTING = 'Listing'

class ResultStatusCode(Enum):
    """
    Python Enum representing the ResultStatusCode type defined in SSC.xsd.
    """
    SUCCESS = 'Success'
    CONDITIONAL_SUCCESS = 'ConditionalSuccess'
    ERROR = 'Error'

class ResultStatusSubCode(Enum):
    """
    Python Enum representing the ResultStatusSubCode type defined in
    SSC.xsd.
    """
    SUCCESS = 'Success'
    MISSING_REQUEST = 'MissingRequest'
    MISSING_SATELLITES = 'MissingSatellites'
    INVALID_BEGIN_TIME = 'InvalidBeginTime'
    INVALID_END_TIME = 'InvalidEndTime'
    INVALID_SATELLITE = 'InvalidSatellite'
    INVALID_TIME_RANGE = 'InvalidTimeRange'
    INVALID_RESOLUTION_FACTOR = 'InvalidResolutionFactor'
    MISSING_OUTPUT_OPTIONS = 'MissingOutputOptions'
    MISSING_COORD_OPTIONS = 'MissingCoordOptions'
    MISSING_COORD_SYSTEM = 'MissingCoordSystem'
    INVALID_COORD_SYSTEM = 'InvalidCoordSystem'
    MISSING_COORD_COMPONENT = 'MissingCoordComponent'
    MISSING_GRAPH_OPTIONS = 'MissingGraphOptions'
    MISSING_COORDINATE_SYSTEM = 'MissingCoordinateSystem'
    MISSING_COORDINATE_COMPONENT = 'MissingCoordinateComponent'
    SERVER_ERROR = 'ServerError'

class SpaceRegion(Enum):
    """
    Python Enum representing the SpaceRegion type defined in SSC.xsd.
    """
    INTERPLANETARY_MEDIUM = 'InterplanetaryMedium'
    DAYSIDE_MAGNETOSHEATH = 'DaysideMagnetosheath'
    NIGHTSIDE_MAGNETOSHEATH = 'NightsideMagnetosheath'
    DAYSIDE_MAGNETOSPHERE = 'DaysideMagnetosphere'
    NIGHTSIDE_MAGNETOSPHERE = 'NightsideMagnetosphere'
    PLASMA_SHEET = 'PlasmaSheet'
    TAIL_LOBE = 'TailLobe'
    LOW_LATITUDE_BOUNDARY_LAYER = 'LowLatitudeBoundaryLayer'
    HIGH_LATITUDE_BOUNDARY_LAYER = 'HighLatitudeBoundaryLayer'
    DAYSIDE_PLASMASPHERE = 'DaysidePlasmasphere'
    NIGHTSIDE_PLASMASPHERE = 'NightsidePlasmasphere'

class SpaceRegionType(Enum):
    """
    Python Enum representing the SpaceRegionType type defined in SSC.xsd.
    """
    INTERPLANETARY_MEDIUM = 'InterplanetaryMedium'
    DAYSIDE_MAGNETOSHEATH = 'DaysideMagnetosheath'
    NIGHTSIDE_MAGNETOSHEATH = 'NightsideMagnetosheath'
    DAYSIDE_MAGNETOSPHERE = 'DaysideMagnetosphere'
    NIGHTSIDE_MAGNETOSPHERE = 'NightsideMagnetosphere'
    PLASMA_SHEET = 'PlasmaSheet'
    TAIL_LOBE = 'TailLobe'
    LOW_LATITUDE_BOUNDARY_LAYER = 'LowLatitudeBoundaryLayer'
    HIGH_LATITUDE_BOUNDARY_LAYER = 'HighLatitudeBoundaryLayer'
    DAYSIDE_PLASMASPHERE = 'DaysidePlasmasphere'
    NIGHTSIDE_PLASMASPHERE = 'NightsidePlasmasphere'

class TimeFormat(Enum):
    """
    Python Enum representing the TimeFormat type defined in SSC.xsd.
    """
    HH_HHHH = 'hh_hhhh'
    HH_MM_SS = 'hh_mm_ss'
    HH_MM = 'hh_mm'

class TraceCoordinateSystem(Enum):
    """
    Python Enum representing the TraceCoordinateSystem type defined
    in SSC.xsd.
    """
    GEO = 'Geo'
    GM = 'Gm'

class TraceType(Enum):
    """
    Python Enum representing the TraceType defined in SSC.xsd.
    """
    B_FIELD = 'BField'
    RADIAL = 'Radial'


#class SscWs:
#    """
#    Class representing the web service interface to NASA's
#    Satelite Situation Center (SSC) <https://sscweb.gsfc.nasa.gov/>.
#    """
#    def __init__(self, endpoint=None, cache=None, timeout=None,
#                 proxy_info=httplib2.proxy_info_from_environment,
#                 ca_certs=None, disable_ssl_certificate_validation=False):
#        self.logger = logging.getLogger(type(self).__name__)
#        self.logger.addHandler(logging.NullHandler())
#
#        self.logger.debug('endpoint = %s', endpoint)
#        self.logger.debug('ca_certs = %s', ca_certs)
#        self.logger.debug('disable_ssl_certificate_validation = %s',
#                          disable_ssl_certificate_validation)
#
#        if endpoint is None:
#            self._endpoint = 'https://sscweb.sci.gsfc.nasa.gov/WS/sscr/2/'
#        else
#            self._endpoint = endpoint
#        self._user_agent = 'sscws/' + __version__ + ' (' + \
#            platform.python_implementation() + ' ' \
#            + platform.python_version() + '; '+ platform.platform() + ')'
#        self._request_headers = {
#            'Content-Type' : 'application/xml',
#            'Accept' : 'application/xml',
#            'User-Agent' : self._user_agent
#        }
#        self._client = httplib2.Http(cache=cache, timeout=timeout,
#                           proxy_info=proxy_info, ca_certs=ca_certs,
#                           disable_ssl_certificate_validation=
#                           disable_ssl_certificate_validation)
#
#    def get_observatories(self):
#        """
#        Gets a description of the available SSC observatories.
#
#        :return: array of dictionaries describing the observatories
#            available at SSC.
#        """
#        return get_observatories()
#
#    def get_locations(self, request):
#        """
#        Gets the given locations DataRequest.
#
#        :param request: dict representation of DataRequest as
#            described in SSC.xsd.
#        :return: dict representation of Result as described in SSC.xsd.
#        """
#        return get_locations(request)



ENDPOINT = "https://sscweb.sci.gsfc.nasa.gov/WS/sscr/2/"
#ENDPOINT = "https://sscweb-dev.sci.gsfc.nasa.gov/WS/sscr/2/"
#ENDPOINT = "http://localhost:8383/WS/sscr/2/"
CLIENT = httplib2.Http(os.path.expanduser('~') + '/.cache',
             disable_ssl_certificate_validation=True)

USER_AGENT = 'sscws/' + __version__ + ' (' + platform.python_implementation() + ' ' \
             + platform.python_version() + '; '+ platform.platform() + ')'

REQUEST_HEADERS = {
    'Content-Type' : 'application/xml',
    'Accept' : 'application/xml',
    'User-Agent' : USER_AGENT
}


def get_observatories():
    """
    Gets a description of the available SSC observatories.

    :return: array of dictionaries describing the observatories
        available at SSC.
    """
#    print(ENDPOINT)
#    print(REQUEST_HEADERS)
    
    dummy_headers, xml = CLIENT.request(ENDPOINT + "observatories", "GET",
                                        headers=REQUEST_HEADERS)
    #print(dummy_headers)
    #print("%s" % xml)

    observatory_response = ET.fromstring(xml)

    observatories = []

    for observatory in observatory_response.findall(\
            '{http://sscweb.gsfc.nasa.gov/schema}Observatory'):

        observatories.append({
            'Id': observatory.find(\
                '{http://sscweb.gsfc.nasa.gov/schema}Id').text,
            'Name': observatory.find(\
                '{http://sscweb.gsfc.nasa.gov/schema}Name').text,
            'Resolution': int(observatory.find(\
                '{http://sscweb.gsfc.nasa.gov/schema}Resolution').text),
            'StartTime': dateutil.parser.parse(observatory.find(\
                '{http://sscweb.gsfc.nasa.gov/schema}StartTime').text),
            'EndTime': dateutil.parser.parse(observatory.find(\
                '{http://sscweb.gsfc.nasa.gov/schema}EndTime').text),
            'ResourceId': observatory.find(\
                '{http://sscweb.gsfc.nasa.gov/schema}ResourceId').text
        })

    return observatories

def sample_orbit(Sat,StartDate,EndDate):
    """
    Obtain 2 days of orbit data and determine the spatial distance range of the orbit
    dict with samplerate, minimum and maximum radial distance from Earth's center.
    """
    data_request = {
        'Description': 'Example locator request.',
        'TimeInterval': {
            'Start': StartDate,
            'End': EndDate
        },
        'Satellites': [],
        'OutputOptions': {
            'AllLocationFilters': True,
            'CoordinateOptions': [],
            'RegionOptions': {
                'Spacecraft': True,
                'RadialTracedFootpoint': False,
                'NorthBTracedFootpoint': False,
                'SouthBTracedFootpoint': False
            },
            'ValueOptions': {
                'RadialDistance': True,
                'BFieldStrength': False,
                'DipoleLValue': False,
                'DipoleInvLat': False
            },
            'DistanceFromOptions': {
                'NeutralSheet': True,
                'BowShock': True,
                'MPause': True,
                'BGseXYZ': True
            },
            'BFieldTraceOptions': [{
                'CoordinateSystem': CoordinateSystem.GEO.value,
                'Hemisphere': Hemisphere.NORTH.value,
                'FootpointLatitude': False,
                'FootpointLongitude': False,
                'FieldLineLength': False
            },{
                'CoordinateSystem': CoordinateSystem.GEO.value,
                'Hemisphere': Hemisphere.SOUTH.value,
                'FootpointLatitude': False,
                'FootpointLongitude': False,
                'FieldLineLength': False,
            }],
        },
        'RegionFilterOptions': {
            'SpaceRegions': {
                'InterplanetaryMedium': True,
                'DaysideMagnetosheath': True,
                'NightsideMagnetosheath': True,
                'DaysideMagnetosphere': True,
                'NightsideMagnetosphere': True,
                'PlasmaSheet': True,
                'TailLobe': True,
                'HighLatitudeBoundaryLayer': True,
                'LowLatitudeBoundaryLayer': True,
                'DaysidePlasmasphere': True,
                'NightsidePlasmasphere': True
            },
            'RadialTraceRegions': {
                'Cusp': {
                    'North': True,
                    'South': True
                },
                'Cleft': {
                    'North': True,
                    'South': True
                },
                'AuroralOval': {
                    'North': True,
                    'South': True
                },
                'PolarCap': {
                    'North': True,
                    'South': True
                },
                'MidLatitude': {
                    'North': True,
                    'South': True
                },
                'LowLatitude': True
            },
            'MagneticTraceRegions': {
                'Cusp': {
                    'North': True,
                    'South': True
                },
                'Cleft': {
                    'North': True,
                    'South': True
                },
                'AuroralOval': {
                    'North': True,
                    'South': True
                },
                'PolarCap': {
                    'North': True,
                    'South': True
                },
                'MidLatitude': {
                    'North': True,
                    'South': True
                },
                'LowLatitude': True
            }            
        }
    }
    data_request['Satellites'].append({
       'Id': Sat,
       'ResolutionFactor': 1
    })
    for component in [CoordinateComponent.X.value,
                      CoordinateComponent.Y.value,
                      CoordinateComponent.Z.value,
                      CoordinateComponent.LAT.value,
                      CoordinateComponent.LON.value,
                      CoordinateComponent.LOCAL_TIME.value]:
        data_request['OutputOptions']['CoordinateOptions'].append({
            'CoordinateSystem': CoordinateSystem.GSE.value,
            'Component': component
        })
    data_request['OutputOptions']['MinMaxPoints'] = 2
#    data_request['OutputOptions']['MinMaxPoints'] = 2

    data_locations=get_locations(data_request)
    return data_locations
                                          
     
def create_example_request():
    """
    Create an example DataRequest.
    >>> currently an incomplete prototype <<<

    :return: dict representation of DataRequest as described in SSC.xsd.
    """
    data_request = {
        'Description': 'Example locator request.',
        'TimeInterval': {
            'Start': datetime(2008, 1, 2, 11, 0, 0, tzinfo=timezone.utc),
            'End': datetime(2008, 1, 2, 11, 59, 59, tzinfo=timezone.utc)
        },
        'BFieldModel': {
            'InternalBFieldModel': InternalBFieldModel.IGRF,
            'ExternalBFieldModel': {
                'Name': ExternalBFieldModel.TSYGANENKO89C,
                'KeyParameterValues': Tsyganenko89cKp.KP_3_3_3
            }
        },
        'Satellites': [],
        'OutputOptions': {
            'AllLocationFilters': True,
            'CoordinateOptions': [],
            'RegionOptions': {
                'Spacecraft': True,
                'RadialTracedFootpoint': True,
                'NorthBTracedFootpoint': True,
                'SouthBTracedFootpoint': True
            },
            'ValueOptions': {
                'RadialDistance': True,
                'BFieldStrength': True,
                'DipoleLValue': True,
                'DipoleInvLat': True
            },
            'DistanceFromOptions': {
                'NeutralSheet': True,
                'BowShock': True,
                'MPause': True,
                'BGseXYZ': True
            },
            'BFieldTraceOptions': [{
                'CoordinateSystem': CoordinateSystem.GEO.value,
                'Hemisphere': Hemisphere.NORTH.value,
                'FootpointLatitude': True,
                'FootpointLongitude': True,
                'FieldLineLength': True
            },{
                'CoordinateSystem': CoordinateSystem.GEO.value,
                'Hemisphere': Hemisphere.SOUTH.value,
                'FootpointLatitude': True,
                'FootpointLongitude': True,
                'FieldLineLength': True,
            }]
        },
        'RegionFilterOptions': {
            'SpaceRegions': {
                'InterplanetaryMedium': True,
                'DaysideMagnetosheath': True,
                'NightsideMagnetosheath': True,
                'DaysideMagnetosphere': True,
                'NightsideMagnetosphere': True,
                'PlasmaSheet': True,
                'TailLobe': True,
                'HighLatitudeBoundaryLayer': True,
                'LowLatitudeBoundaryLayer': True,
                'DaysidePlasmasphere': True,
                'NightsidePlasmasphere': True
            },
            'RadialTraceRegions': {
                'Cusp': {
                    'North': True,
                    'South': True
                },
                'Cleft': {
                    'North': True,
                    'South': True
                },
                'AuroralOval': {
                    'North': True,
                    'South': True
                },
                'PolarCap': {
                    'North': True,
                    'South': True
                },
                'MidLatitude': {
                    'North': True,
                    'South': True
                },
                'LowLatitude': True
            },
            'MagneticTraceRegions': {
                'Cusp': {
                    'North': True,
                    'South': True
                },
                'Cleft': {
                    'North': True,
                    'South': True
                },
                'AuroralOval': {
                    'North': True,
                    'South': True
                },
                'PolarCap': {
                    'North': True,
                    'South': True
                },
                'MidLatitude': {
                    'North': True,
                    'South': True
                },
                'LowLatitude': True
            }
        }
    }
    for sat in ['themisa', 'spase://SMWG/Observatory/THEMIS/B']:
        data_request['Satellites'].append({
            'Id': sat,
            'ResolutionFactor': 2
        })
    for component in [CoordinateComponent.X.value,
                      CoordinateComponent.Y.value,
                      CoordinateComponent.Z.value,
                      CoordinateComponent.LAT.value,
                      CoordinateComponent.LON.value,
                      CoordinateComponent.LOCAL_TIME.value]:
        data_request['OutputOptions']['CoordinateOptions'].append({
            'CoordinateSystem': CoordinateSystem.GSE.value,
            'Component': component
        })

    data_request['OutputOptions']['MinMaxPoints'] = 2

    data_locations=get_locations(data_request)
                                      
    return data_locations


# pylint: disable=too-many-statements

def __create_xml_data_request(data_request):
    """
    Creates an XML DataRequest XML string from a dict representation
    of a DataRequest.
    >>> currently an incomplete prototype <<<

    :return: string containing the XML representation of the given
        dict representation of a DataRequest as described in SSC.xsd.
    """

    builder = ET.TreeBuilder()
    builder.start('DataRequest', {'xmlns': 'http://sscweb.gsfc.nasa.gov/schema'})
    builder.start('Description', {})
    builder.data(data_request['Description'])
    builder.end('Description')
    builder.start('TimeInterval', {})
    builder.start('Start', {})
    builder.data(data_request['TimeInterval']['Start'].isoformat())
    builder.end('Start')
    builder.start('End', {})
    builder.data(data_request['TimeInterval']['End'].isoformat())
    builder.end('End')
    builder.end('TimeInterval')

    for sat in data_request['Satellites']:
        builder.start('Satellites', {})
        builder.start('Id', {})
        builder.data(sat['Id'])
        builder.end('Id')
        builder.start('ResolutionFactor', {})
        builder.data(str(sat['ResolutionFactor']))
        builder.end('ResolutionFactor')
        builder.end('Satellites')

    if 'BFieldModel' in data_request:
        builder.start('BFieldModel', {})
        builder.start('InternalBFieldModel', {})
        builder.data(data_request['BFieldModel']['InternalBFieldModel'].value)
        builder.end('InternalBFieldModel')
        external_model = data_request['BFieldModel']['ExternalBFieldModel']
        external_model_name = external_model['Name'].value
        builder.start('ExternalBFieldModel', {
            'xmlns:xsi': 'http://www.w3.org/2001/XMLSchema-instance',
            'xsi:type': external_model_name
        })
        if external_model_name is ExternalBFieldModel.TSYGANENKO96:
            builder.start('SolarWindPressure', {})
            builder.data(str(external_model['SolarWindPressure']))
            builder.end('SolarWindPressure')
            builder.start('DstIndex', {})
            builder.data(str(external_model['DstIndex']))
            builder.end('DstIndex')
            builder.start('ByImf', {})
            builder.data(str(external_model['ByImf']))
            builder.end('ByImf')
            builder.start('BzImf', {})
            builder.data(str(external_model['BzImf']))
            builder.end('BzImf')
        else:
            builder.start('KeyParameterValues', {})
            builder.data(external_model['KeyParameterValues'].value)
            builder.end('KeyParameterValues')
        builder.end('ExternalBFieldModel')
        builder.end('BFieldModel')

    builder.start('OutputOptions', {})
    builder.start('AllLocationFilters', {})
    builder.data(str(data_request['OutputOptions']['AllLocationFilters']).lower())
    builder.end('AllLocationFilters')

    for coord_option in data_request['OutputOptions']['CoordinateOptions']:
        builder.start('CoordinateOptions', {})
        builder.start('CoordinateSystem', {})
        builder.data(coord_option['CoordinateSystem'])
        builder.end('CoordinateSystem')
        builder.start('Component', {})
        builder.data(coord_option['Component'])
        builder.end('Component')
        builder.end('CoordinateOptions')

    if 'RegionOptions' in data_request['OutputOptions']:
        region_option = data_request['OutputOptions']['RegionOptions']
        builder.start('RegionOptions', {})
        builder.start('Spacecraft', {})
        builder.data(str(region_option['Spacecraft']).lower())
        builder.end('Spacecraft')
        builder.start('RadialTracedFootpoint', {})
        builder.data(str(region_option['RadialTracedFootpoint']).lower())
        builder.end('RadialTracedFootpoint')
        builder.start('NorthBTracedFootpoint', {})
        builder.data(str(region_option['NorthBTracedFootpoint']).lower())
        builder.end('NorthBTracedFootpoint')
        builder.start('SouthBTracedFootpoint', {})
        builder.data(str(region_option['SouthBTracedFootpoint']).lower())
        builder.end('SouthBTracedFootpoint')
        builder.end('RegionOptions')

    if 'ValueOptions' in data_request['OutputOptions']:
        value_option = data_request['OutputOptions']['ValueOptions']
        builder.start('ValueOptions', {})
        builder.start('RadialDistance', {})
        builder.data(str(value_option['RadialDistance']).lower())
        builder.end('RadialDistance')
        builder.start('BFieldStrength', {})
        builder.data(str(value_option['BFieldStrength']).lower())
        builder.end('BFieldStrength')
        builder.start('DipoleLValue', {})
        builder.data(str(value_option['DipoleLValue']).lower())
        builder.end('DipoleLValue')
        builder.start('DipoleInvLat', {})
        builder.data(str(value_option['DipoleInvLat']).lower())
        builder.end('DipoleInvLat')
        builder.end('ValueOptions')

    if 'DistanceFromOptions' in data_request['OutputOptions']:
        distance_option = data_request['OutputOptions']['DistanceFromOptions']
        builder.start('DistanceFromOptions', {})
        builder.start('NeutralSheet', {})
        builder.data(str(distance_option['NeutralSheet']).lower())
        builder.end('NeutralSheet')
        builder.start('BowShock', {})
        builder.data(str(distance_option['BowShock']).lower())
        builder.end('BowShock')
        builder.start('MPause', {})
        builder.data(str(distance_option['MPause']).lower())
        builder.end('MPause')
        builder.start('BGseXYZ', {})
        builder.data(str(distance_option['BGseXYZ']).lower())
        builder.end('BGseXYZ')
        builder.end('DistanceFromOptions')

    for bfield_trace_option in data_request['OutputOptions']['BFieldTraceOptions']:
        builder.start('BFieldTraceOptions', {})
        builder.start('CoordinateSystem', {})
        builder.data(bfield_trace_option['CoordinateSystem'])
        builder.end('CoordinateSystem')
        builder.start('Hemisphere', {})
        builder.data(bfield_trace_option['Hemisphere'])
        builder.end('Hemisphere')
        builder.start('FootpointLatitude', {})
        builder.data(str(bfield_trace_option['FootpointLatitude']).lower())
        builder.end('FootpointLatitude')
        builder.start('FootpointLongitude', {})
        builder.data(str(bfield_trace_option['FootpointLongitude']).lower())
        builder.end('FootpointLongitude')
        builder.start('FieldLineLength', {})
        builder.data(str(bfield_trace_option['FieldLineLength']).lower())
        builder.end('FieldLineLength')
        builder.end('BFieldTraceOptions')

    builder.start('MinMaxPoints', {})
    builder.data(str(data_request['OutputOptions']['MinMaxPoints']))
    builder.end('MinMaxPoints')
    builder.end('OutputOptions')
    builder.end('DataRequest')

    xml_data_request = builder.close()

    return xml_data_request

# pylint: enable=too-many-statements



def get_locations(data_request):
    """
    Gets the given locations DataRequest.

    :param request: dict representation of DataRequest as
        described in SSC.xsd.
    :return: dict representation of Result as described in SSC.xsd.
    """

    xml_data_request = __create_xml_data_request(data_request)

    #print(ET.tostring(xml_data_request))

    dummy_headers, xml = CLIENT.request(ENDPOINT + "locations", "POST",
                                        body=ET.tostring(xml_data_request),
                                        headers=REQUEST_HEADERS)
    #print("%s" % xml)

    result_element = ET.fromstring(xml).find(\
                         '{http://sscweb.gsfc.nasa.gov/schema}Result')

    return __get_result(result_element)


def __get_result(result_element):
    """
    Creates a dict representation of a DataResult from an ElementTree
    representation
    >>> currently an incomplete prototype <<<

    :return: dict representation of the given ElementTree DataResult
        as described in SSC.xsd.
    """

    result = {
        'StatusCode': result_element.find(\
               '{http://sscweb.gsfc.nasa.gov/schema}StatusCode').text,
        'StatusSubCode': result_element.find(\
               '{http://sscweb.gsfc.nasa.gov/schema}StatusSubCode').text,
        'Data': []
    }

    data_i = -1

    for data_element in result_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}Data'):

        data_i += 1

        coords_element = data_element.find(\
               '{http://sscweb.gsfc.nasa.gov/schema}Coordinates')
        coordinates = {
            'CoordinateSystem': CoordinateSystem(coords_element.find(\
               '{http://sscweb.gsfc.nasa.gov/schema}CoordinateSystem').text),
            'X': [],
            'Y': [],
            'Z': [],
            'Latitude': [],
            'Longitude': [],
            'LocalTime': []
        }
        result['Data'].append({
            'Id': data_element.find(\
                      '{http://sscweb.gsfc.nasa.gov/schema}Id').text,
            'Coordinates': coordinates,
            'Time': [],
            'BTraceData': [],
            'RadialLength': [],
            'MagneticStrength': [],
            'NeutralSheetDistance': [],
            'BowShockDistance': [],
            'MagnetoPauseDistance': [],
            'DipoleLValue': [],
            'DipoleInvariantLatitude': [],
            'SpacecraftRegion': [],
            'RadialTracedFootpointRegions': [],
            'BGseX': [],
            'BGseY': [],
            'BGseZ': [],
            'NorthBTracedFootpointRegions': [],
            'SouthBTracedFootpointRegions': []
        })

        for x_coord in coords_element.findall(\
                '{http://sscweb.gsfc.nasa.gov/schema}X'):

            result['Data'][data_i]['Coordinates']['X'].append(\
                float(x_coord.text))

        for y_coord in coords_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}Y'):

            result['Data'][data_i]['Coordinates']['Y'].append(\
                float(y_coord.text))

        for z_coord in coords_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}Z'):

            result['Data'][data_i]['Coordinates']['Z'].append(\
                float(z_coord.text))

        for lat_coord in coords_element.findall(\
                '{http://sscweb.gsfc.nasa.gov/schema}Latitude'):

            result['Data'][data_i]['Coordinates']['Latitude'].append(\
                float(lat_coord.text))

        for lon_coord in coords_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}Longitude'):

            result['Data'][data_i]['Coordinates']['Longitude'].append(\
                float(lon_coord.text))

        for lt_coord in coords_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}LocalTime'):

            result['Data'][data_i]['Coordinates']['LocalTime'].append(\
                float(lt_coord.text))

        for time in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}Time'):

            result['Data'][data_i]['Time'].append(\
                dateutil.parser.parse(time.text))

        for b_trace_data in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}BTraceData'):

             result['Data'][data_i]['BTraceData'].append({
                 'CoordinateSystem': CoordinateSystem(b_trace_data.find(\
                     '{http://sscweb.gsfc.nasa.gov/schema}CoordinateSystem').text),
                 'Hemisphere': Hemisphere(b_trace_data.find(\
                     '{http://sscweb.gsfc.nasa.gov/schema}Hemisphere').text),
                 'Latitude': [], 
                 'Longitude': [],
                 'ArcLength': []
             })
             for lat in b_trace_data.findall(\
                 '{http://sscweb.gsfc.nasa.gov/schema}Latitude'):

                 result['Data'][data_i]['BTraceData'][-1]['Latitude'].append(\
                     float(lat.text))

             for lon in b_trace_data.findall(\
                 '{http://sscweb.gsfc.nasa.gov/schema}Longitude'):

                 result['Data'][data_i]['BTraceData'][-1]['Longitude'].append(\
                     float(lon.text))

             for arc_length in b_trace_data.findall(\
                 '{http://sscweb.gsfc.nasa.gov/schema}ArcLength'):

                 result['Data'][data_i]['BTraceData'][-1]['ArcLength'].append(\
                     float(arc_length.text))

        for radial_length in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}RadialLength'):

            result['Data'][data_i]['RadialLength'].append(\
                float(radial_length.text))

        for magnetic_strength in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}MagneticStrength'):

            result['Data'][data_i]['MagneticStrength'].append(\
                float(magnetic_strength.text))

        for neutral_sheet_distance in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}NeutralSheetDistance'):

            result['Data'][data_i]['NeutralSheetDistance'].append(\
                float(neutral_sheet_distance.text))

        for bow_shock_distance in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}BowShockDistance'):

            result['Data'][data_i]['BowShockDistance'].append(\
                float(bow_shock_distance.text))

        for magneto_pause_distance in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}MagnetoPauseDistance'):

            result['Data'][data_i]['MagnetoPauseDistance'].append(\
                float(magneto_pause_distance.text))

        for dipole_l_value in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}DipoleLValue'):

            result['Data'][data_i]['DipoleLValue'].append(\
                float(dipole_l_value.text))

        for dipole_invariant_latitude in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}DipoleInvariantLatitude'):

            result['Data'][data_i]['DipoleInvariantLatitude'].append(\
                float(dipole_invariant_latitude.text))

        for spacecraft_region in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}SpacecraftRegion'):

            result['Data'][data_i]['SpacecraftRegion'].append(\
                SpaceRegion(spacecraft_region.text))

        for radial_footpoint_region in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}RadialTracedFootpointRegions'):

            result['Data'][data_i]['RadialTracedFootpointRegions'].append(\
                FootpointRegion(radial_footpoint_region.text))

        for b_gse_x in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}BGseX'):

            result['Data'][data_i]['BGseX'].append(\
                float(b_gse_x.text))

        for b_gse_y in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}BGseY'):

            result['Data'][data_i]['BGseY'].append(\
                float(b_gse_y.text))

        for b_gse_z in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}BGseZ'):

            result['Data'][data_i]['BGseZ'].append(\
                float(b_gse_z.text))

        for b_traced_footpoint_region in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}NorthBTracedFootpointRegions'):

            result['Data'][data_i]['NorthBTracedFootpointRegions'].append(\
                FootpointRegion(b_traced_footpoint_region.text))

        for b_traced_footpoint_region in data_element.findall(\
               '{http://sscweb.gsfc.nasa.gov/schema}SouthBTracedFootpointRegions'):

            result['Data'][data_i]['SouthBTracedFootpointRegions'].append(\
                FootpointRegion(b_traced_footpoint_region.text))

    #print(result)

    return result


def print_locations_result(result):
    """
    Prints a Result document.

    :param result: dict representation of Result as described in SSC.xsd.
    """

    print('StatusCode:', result['StatusCode'],
          'StatusSubCode:', result['StatusSubCode'])

    for data in result['Data']:
        coords = data['Coordinates']
        print(data['Id'], coords['CoordinateSystem'].value)
        print('Time                     ', 'X                     ',
              'Y                     ', 'Z                     ')
        for index in range(len(data['Time'])):
            print(data['Time'][index], coords['X'][index],
                  coords['Y'][index], coords['Z'][index])

        for b_trace in data['BTraceData']:

            print(b_trace['CoordinateSystem'].value, 
                  b_trace['Hemisphere'].value)
            print('Time                          ', 'Latitude        ',
                  'Longitude   ', 'Arc Length')
            for index in range(len(data['Time'])):
                print(data['Time'][index], 
                      '{:15.5f} {:15.5f} {:15.5f}'.format(\
                          b_trace['Latitude'][index],
                          b_trace['Longitude'][index], 
                          b_trace['ArcLength'][index]))

        if 'RadialLength' in data:
            print('Time                     ', 'Radial Length         ')
            for index in range(len(data['Time'])):
                print(data['Time'][index], data['RadialLength'][index])

        # repeat for MagneticStrength, NeutralSheetDistance, 
        # BowShockDistance, MagnetoPauseDistance, DipoleLValue,
        # DipoleInvariantLatitude, SpacecraftRegion, 
        # RadialTracedFootpointRegions

        if 'BGseX' in data:

            print('Time                     ', 'B Strength GSE        ')
            print('                             X                     ', 
                  'Y                     ', 'Z')
            for index in range(len(data['Time'])):
                print(data['Time'][index], data['BGseX'][index],
                      data['BGseY'][index], 
                      data['BGseZ'][index])

        if 'NorthBTracedFootpointRegion' in data and \
           'SouthBTracedFootpointRegion' in data:

            print('                 B-Traced Footpoint Region')
            print('Time                     ', 'North            ', 
                  'South           ')
            for index in range(len(data['Time'])):
                print(data['Time'][index], 
                      data['NorthBTracedFootpointRegion'][index].value,
                      data['SouthBTracedFootpointRegion'][index].value)

def get_orbit_regions(orbit_data):
# calculate radial distances, minimum, maximum and regions encountered
# radial distances from Earth's center defining regions (disciplines)
# numerical values taken from SSCWEB_satellite_inventory.pl script that we are replacing
    rmax_it=1.2  # ionosphere-thermosphere
    rmin_im=1.2  # inner magnetosphere
    rmax_im=10.
    rmin_gm=3.   # magnetosphere
    rmax_gm=65   
    rmin_ih=30   # (inner) heliosphere

    X=np.array(orbit_data['Data'][0]['Coordinates']['X']) #
    Y=np.array(orbit_data['Data'][0]['Coordinates']['Y']) # positions im km in GSE
    Z=np.array(orbit_data['Data'][0]['Coordinates']['Z']) #
    RE=6371.200 # R_E in km
    R2=(X*X+Y*Y+Z*Z)/(RE*RE) # square of radial distance from Earth's center [R_E^2]
#    R=np.array(orbit_data['Data'][0]['Coordinates']['RadialLength'])/RE
    rmin=math.sqrt(R2.min())
    rmax=math.sqrt(R2.max())
    sat_discipline=""
    if rmin <= rmax_it:
        sat_discipline="%s,IT" % sat_discipline
    if rmax >= rmin_im and rmin < rmax_im:
        sat_discipline="%s,IM" % sat_discipline
    if rmax >= rmin_gm and rmin < rmax_gm:
        sat_discipline="%s,GM" % sat_discipline
    if rmax >= rmin_ih:
        sat_discipline="%s,SH" % sat_discipline

    sat_discipline=sat_discipline[1:len(sat_discipline)];
# sometimes one only gets a single data point (IMAGE satellite)
    return {'rmin': rmin, 'rmax': rmax, 'regions': sat_discipline}

# deleted main() routine
