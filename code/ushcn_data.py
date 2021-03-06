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
which contains meta-data ab
out the Station from where the raw data comes.

In the future, an additional class will be added to contain collections of
'Series' data objects.

"""
__docformat__ = "restructuredtext"

from util import compute_monthly_anomalies

#: The element codes identifying what variable is contained in a USHCN dataset,
#: and the scaling factor to convert the recorded values to observations.
ELEM_CODES = { 1: { 'type': 'max', 'scale': 0.1 },
               2: { 'type': 'min', 'scale': 0.1 },
               3: { 'type': 'avg', 'scale': 0.1 },
               4: { 'type': 'pcp', 'scale': 0.01 } }

#: Possible types of USHCN data.
ELEM_TYPES = [data['type'] for data in ELEM_CODES.values()]

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
        return "%s (coop_id=%6s)" % (self.name, self.coop_id)
    
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
        
        MISSING_VAL = MISSING
        variable = None
        variable_str = None
        if 'variable' in k:
            variable = k['variable']
            
            variable_info = ELEM_CODES[variable]
            MISSING_VAL = MISSING
            variable_str = variable_info['type']            
                
            del k['variable']
            
            self.MISSING_VAL = MISSING_VAL
            self._variable = variable
            self._variable_str = variable_str

        series = None
        years = None
        subset_years = None
        self._deleted_months = []
        if 'series' in k:
            if not 'years' in k:
                raise MissingDataError("years")
            series = k['series']
            years = k['years']
            if 'subset_years' in k:
                subset_years = k['subset_years']
            
            ## We have data, but there's probably missing years in there. We
            ## need to go ahead and fill in missing data so that we have 
            ## semi-continuous data (or at least, continuous placeholders for
            ## data) running from years[0] to years[-1].
            series, years = self._fill_missing(series, years, self.MISSING_VAL,
                                               subset_years)
            
            del k['series']
            del k['years']
            self.set_series(series, years)
                    
        self.__dict__.update(k)        
            
    def __repr__(self):
        # Assume that this series is associated with a station and knows that
        # station's 6-digit USHCN Coop Id and name
        return "%s (coop_id=%6s)" % (self.name, self.coop_id)
    
    @property
    def variable(self):
        """Get the string indicating what type of data is contained in this
        series."""
        return self._variable_str
    
    @property
    def series(self):
        """Get the actual data contained in this series."""
        return self._series
    
    @property 
    def deleted_months(self):
        """Get the months which have been tagged for deleted"""
        return self._deleted_months
        
    @property 
    def __len__(self):
        """The length of the series data contained in this object."""
        return len(self._series)
    
    @property
    def first_year(self):
        """The year in which the data in this series begins."""
        return self._years[0]
    
    @property
    def last_year(self):
        """The last year for which there is data in this series."""
        return self._years[-1]
        
    def trunc_series(self, begyr, endyr):
        trunc_series = []
#        new_years = range(begyr, endyr)
#        for year in new_years:
#            if year in self._years:
#                year_idx = self._years.index(year)
#                trunc_series.append(self._data[year_idx])
#            else:
#                trunc_series.append()
        for (year, data) in zip(self._years, self._series):
            if (begyr <= year and year < endyr):
                trunc_series.append(data)
        return trunc_series
        
    def set_series(self, series, years):
        """Set the actual data series in this object."""
        self._years = list(years)
        self._series = list(series)
        self._monthly = self._flatten_months(series)
        
    def delete_months(self, months):
        """Add the given months to the list of deleted months"""
        self._deleted_months.extend(months)
        
    @property
    def years(self):
        """Returns the list of years corresponding to the data in this Series
        instance"""
        return self._years        
        
    @property
    def monthly_series(self):
        """Returns a flattened, monthly list of the data values in this object,
        of length len(years)*12 as opposed to len(years)."""
        return self._monthly
    
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
    
    def _fill_missing(self, series, years, fill_val=-9999, subset_years=None):
        """Determines where there are missing years in the provided series
        of data and fills them with the provided fill value.
        """
        missing_yearly = [fill_val]*13
        filled_series = []
        
        if subset_years:
            filled_years = subset_years
        else:
            filled_years = range(years[0], years[-1]+1)
            
        for year in filled_years:
            if year in years:
                filled_series.append(series[years.index(year)])
            else:
                filled_series.append(missing_yearly)
        
        return filled_series, filled_years
    
    @property
    def monthly_anomaly_series(self):
        """Returns the monthly anomalies computed from this series as a flat
        list of length (len(self.years)*12)
        """
        anomalies = compute_monthly_anomalies(self.series, 
                                              self.MISSING_VAL)
        flat_anomalies = self._flatten_months(anomalies)
        
        return flat_anomalies        
    
class Network(object):
    """USHCN Network of stations/associated data.
    
    A 'Network' instance is an object which holds meta-data for a specific
    group of stations from the USHCN and their data, which a user would like
    to homogenize. 
    
    :Ivar stations:
        A dictionary of 'Station' objects, organized as follows:
        stations = dict(station_a.coop_id=station_a,
                        station_b.coop_id=station_b, ... )
    :Ivar raw_series:
        A dictionary of 'Series' objects, organized as follows:
        raw_series = dict(series_a.coop_id=series_a,
                          series_b.coop_id=series_b, ... )
        
    Both series and stations need to have the same set of keys; this will be
    explicitly checked during construction.
    
    ...
    
    :Raises ValueError:
        If stations and raw_series supplied to the constructor have a different
        set of keys.
                    
    """
    
    def __init__(self, stations, raw_series, **k):
        
        # Check if stations and raw_series have the same set of keys
        stations_ids = sorted(stations.keys())
        series_ids = sorted(raw_series.keys())
        if not (stations_ids == series_ids):
            raise ValueError("'stations_ids' and 'series_ids' must have same keys!")
                
        self._stations = stations
        self._raw_series = raw_series
        
        
        self.__dict__.update(k)
    
    def __repr__(self):
        
        out_str = ""
        if hasattr(self, "name"):
            out_str = out_str + ("%s:\n    " % self.name)
            
        stations_list = self._stations.itervalues()
        stations_str = ",\n    ".join(map(str, stations_list))
        
        return out_str + stations_str
        
    @property
    def stations(self):
        return self._stations
    
    @property
    def raw_series(self):
        return self._raw_series
    
    #######################
    
    def add_station(self, station):
        self._series.append(station)
        
    def add_series(self, series):
        self._series.append(series)
        
    def update_series(self, new_series):
        self._raw_series = new_series

class MissingDataError(Exception):
    """Exception raised if a user tries to instantiate a Station or Series 
    object but fails to supply a necessary data field.
    
    :Ivar field:
        The field which was not supplied.
    
    """
    def __init__(self, field):
        self.field = field
        
    def __repr__(self):
        return "Need to supply data for the field %s." % self.field