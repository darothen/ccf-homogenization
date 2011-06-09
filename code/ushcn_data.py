#!/usr/bin/env python
#
# Daniel Rothenberg, 2011-06-06

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
    
    def update_values(self, **values):
        self.__dict__.update(values)
    
    @property
    def values_dict(self):
        return self.__dict__
    
# TODO: This is a work in progress.
class Series(object):
    """Monthly data Series.
    
    An instance of 'Series' contains a series of monthly data (either average,
    minimum, maximum temperatures or precipitation), accessible via the 'series'
    property. 
    
    The following are meta-data which should be expected to accompany a 'Series'
    object:
    
    :Ivar coop_id:
        A 6-digit identification number for the series associated with
        the station which produced it. Should be provided as a string so that
        if the station begins with a '0', it can still be correctly identified
        
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
    
    :Ivar coop_1, coop_2, coop_3:
        The 6 digit Coop Id's for the first[, second, and third] stations in
        chronological order (if applicable) whose records were joined to form
        the longer time series contained here.
    
    :Ivar utc_offset:
        An integer in (-12, 12) giving the time difference between Coordinated
        Universal Time (UTC) and the local standard time at the station
        generating this data.
    
    """
    def __init__(self, **k):
        first_year = None
        if 'first_year' in k:
            first_year = k['first_year']
            del k['first_year']
            self._first_year = first_year

        series = None
        if 'series' in k:
            series = k['series']
            del k['series']
            self.set_series(series, first_year)
        
        MISSING_VAL = MISSING
        variable = None
        if 'variable' in k:
            variable = k['variable']
            if variable in (1, 2, 3):
                MISSING_VAL = MISSING*0.1
            if variable in (4, ):
                MISSING_VAL = MISSING*0.01
            del k['variable']
            self.MISSING_VAL = MISSING_VAL
            self._variable = variable
            
        self.__dict__.update(k)
        
            
    def __repr__(self):
        # Assume that this series is associated with a station and knows that
        # station's 6-digit USHCN Coop Id and name
        return "%s (coop_id=%6s)" % (self.name, self.coop_id)
    
    @property
    def series(self):
        """Get the actual data contained in this series."""
        return self._series
    
    @property 
    def __len__(self):
        """The length of the series data contained in this object."""
        return len(self._series)
    
    @property
    def first_year(self):
        """The year in which the data in this series begins."""
        return self._first_year
    
    @property
    def last_year(self):
        """The last year for which there is data in this series."""
        return (self._first_year + len(self._series) - 1)
        
    def set_series(self, series, first_year):
        """Set the actual data series in this object."""
        self._first_year = first_year
        self._series = list(series)
        
    @property
    def monthly_series(self):
        """Returns a flattened, monthly list of the data values in this object,
        of length len(years)*12 as opposed to len(years)."""
        return self._flatten_months(self._series)
    
    def _flatten_months(self, l):
        """Flattens and extracts a timeseries of monthly datavalues from a
        Series."""
        new_list = []
        for obj in l:
            if isinstance(obj, (tuple, list)):
                if len(obj) == 13:
                    new_list.extend(obj[:-1])
                else:
                    new_list.extend(obj)
            else:
                new_list.append(obj)
        return new_list
                
