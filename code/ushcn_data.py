#!/usr/bin/env python
# $URL$
# $Rev$

"""Classes for USHCN data.

The United States Historical Climatology Network (USHCN, 
http://cdiac.ornl.gov/epubs/ndp/ushcn/access.html) records monthly observations
of maximum, minimum, and average temperatures, and precipitation. The classes
defined here can be used for working with monthly data for any of these
variables. 

The 'Series' class here is the primary tool for holding raw data for any
of the USHCN variables. Any 'Series' object must reference a 'Station' object,
which contains meta-data about the Station from where the raw data comes.

In the future, an additional class will be added to contain collections of
'Series' data objects.

"""
__docformat__ = "restructuredtext"

#: The value used in the USHCN datasets to indicate a missing data value
MISSING = -9999

def invalid(v):
    return v == MISSING

def valid(v):
    return not invalid(v)

class Station(object):
    """A Station object containing station meta-data and raw observations.
    
    This holds meta-data and information about a single USHCN monitoring station.
    Not all of the attributes and data available will be used by the homogenization
    algorithms in this library.
    """
    def __init__(self, **values):
        self.__dict__.update(values)
        
    def __repr__(self): 
        return "Station(%r)" % self.__dict__
    
# TODO: This is a work in progress.
class Series(object):
    """Monthly data Series.
    
    An instance of 'Series' contains a series of monthly data (either average,
    minimum, maximum temperatures or precipitation), accessible via the 'series'
    property. This property should **always** be treated as read-only. Also, 
    other series meta-data is provided by read-only properties. Generally, only
    the defined mutator methods should be used to modify properties.
    
    The following are meta-data which should be expected to accompany a 'Series'
    object:
    
    :Ivar id:
        A 6-digit integer identification number for the series associated with
        the station which produced it.
        
    :Ivar type:
        A 3-letter string ('avg', 'max', 'min', 'pcp') indicating what type of 
        data is contained here.
    
    :Ivar lat, lon:
        Coordinates where the station generating this data is located, in decimal
        degrees.
        
    :Ivar elev:
        The elevation of the station which generated this data, in meters. If 
        the station elevation is missing or unknown, elev = -999.9
    
    :Ivar state:
        The U.S. postal code for the state in which the station generating this 
        data is located.
        
    :Ivar name:
        The name of the station generating this data.
    
    :Ivar comp_1, comp_2, comp_3:
        The 6 digit Coop Id's for the first[, second, and third] stations in
        chronological order (if applicable) whose records were joined to form
        the longer time series contained here.
    
    :Ivar utc_offset:
        An integer in (-12, 12) giving the time difference between Coordinated
        Universal Time (UTC) and the local standard time at the station
        generating this data.
    
    """
    