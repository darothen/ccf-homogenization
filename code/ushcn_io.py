#!/usr/bin/env python
#
# Daniel Rothenberg, 2011-06-08

"""Input/output methods for USHCN datasets.

The data used to construct the United States Historical Climatology Network
(USHCN) is available to freely download from the National Climatic Data Center
at ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/v2/monthly/. The methods contained 
here are provided to help download data for use in the homogenization routines
in this library.

Furthermore, these methods are intended to work with the classes in the
'ushcn_data.py' file, which serve as blueprints for storing and accessing
data once it has been read into this program.

"""
__docformat__ = "restructuredtext"

import os
import urllib
import gzip

from ushcn_data import Station, Series

#: The root directory for the USHCN monthly data ftp site
USHCN_FTP_ROOT = "http://cdiac.ornl.gov/ftp/ushcn_v2_monthly/"

#: The name of the USHCN station list file
USHCN_STATION_LIST = "ushcn-stations.txt"

#: A pattern to use for matching and downloading data from the USHCN_FTP_ROOT
#: site.
USHCN_RAW_DATA_PATTERN = "9641C_200912_%s.%s"

#: The element codes identifying what variable is contained in a USHCN dataset
ELEM_CODES = { 'max':1, 'min':2, 'avg':3, 'pcp':4 }

#: Possible sources of USHCN data.
SOURCES = ['raw', 'tob', 'F52']

#: Possible types of USHCN data.
ELEM_TYPES = ELEM_CODES.keys()

#: Possible flags on USHCN data values.
USHCN_FLAGS = ['E', 'I', 'Q', 'X']

def is_flagged(value):
    """Determines if a raw value carries a USHCN flag, signifying that it was
    estimated in some way."""
    for flag in USHCN_FLAGS:
        if value.endswith(flag):
            return True
    return False

def get_ushcn_data(source, variable, stations=None):    
    """Download USHCN datasets for use in a homogenization algorithm.
    
    :Param source:
        The source dataset provided by USHCN which the user wishes to download;
        must be either "F52" (bias-adjusted mean monthly values with
        estimates for missing values), "tob" (mean monthly values adjusted only
        for the time of observation bias), or "raw" (unadjusted mean monthly
        values).
    :Param variable:
        The variable whose mean monthly values the user wishes to download; must
        be either "max" (monthly maximum temperatures), "min" (monthly minimum
        temperatures), "avg" (average of monthly maximum and minimum
        temperatures), or "pcp" (total monthly precipitation).
    :Param stations:
        (optional) A list stations whose data should be returned. If omitted,
        this method will return all the data it finds.
        
    :Return:
        A collection of USHCN data series.
        
    """
    if source not in SOURCES:
        raise ArgumentError("source", source, proper=SOURCES)
    if variable not in ELEM_TYPES:
        raise ArgumentError("variable", variable, proper=ELEM_TYPES)    
    
    data_fn = (USHCN_RAW_DATA_PATTERN % (source, variable)) + ".gz"
    station_list_fn = USHCN_STATION_LIST
    
    data_dir = os.path.join(os.getcwd(), "data")
    
    def ftp_path(fn):
        return USHCN_FTP_ROOT+fn
    
    def data_path(fn):
        return os.path.join(data_dir, fn)
    
    ## Check if there is a local copy of the USHCN source data file and
    ## station list; if not, download and uncompress it.
    for fn in [data_fn, station_list_fn]:
        if not os.path.exists(data_path(fn)):
            print "Attempting to download %s..." % ftp_path(fn),
            (fname, headers) = urllib.urlretrieve(ftp_path(fn), data_path(fn))
            print "done. Local file is %s" % fname
        else:
            print "Found data/%s" % fn  
            fname = data_path(fn)

        if fname.endswith(".gz"):
            uncompressed_file = gzip.open(fname, 'rb')
            f_out = open(data_path(fname[:-3]), 'wb') 
            f_out.writelines(uncompressed_file.read())
            f_out.close()
    
    station_list_file = None
    try:
        station_list_file = open(data_path(station_list_fn), 'rb')
    except IOError:
        print "Could not open file %s." % (data_path(station_list_fn[:-3]))
        
    all_stations = dict()
    for station_str in station_list_file:
        new_station = read_station_string(station_str)
        all_stations[new_station.coop_id] = new_station
    station_list_file.close()
    
    ## Now that we have all the stations and their corresponding metadata,
    ## we can process the data series associated with the station.    
    if stations:
        stations = sorted(stations)
    else:
        stations = sorted(all_stations.keys())
    
    try:
        big_datafile = open(data_path(data_fn)[:-3], 'rb')
    except IOError:
        print "Could not open file %s." % (data_path(data_fn)[:-3])
    
    appended_stations = []
    for data_str in big_datafile:
        line_info = read_ushcn_dataset_string(data_str)
        (station_id, variable, year, monthly, annual) = line_info
        yearly_data = monthly + [annual]
        
        station = all_stations[station_id]
        
        if station_id in appended_stations:
            station.series = station.series + [yearly_data]
            #all_stations[station_id] = station
        else:
            if station_id in stations:
                station.update_values(series=[yearly_data], first_year=year,
                                      variable=variable)
                appended_stations.append(station_id)
                
    all_series = dict()
    for station_id in stations:
        station = all_stations[station_id]
        series = Series(**station.values_dict)
        all_series[station_id] = series
    return all_series

def read_ushcn_dataset_string(data_str):
    """
    
    """
    data_str = data_str.strip()
    data_bits = data_str.split()
    
    identifier = data_bits[0]
    id = identifier[:6]
    variable = identifier[6]
    year = identifier[7:11]
    
    monthly_values = data_bits[1:13]
    monthly = []
    for val in monthly_values:
        if is_flagged(val):
            stripped_val = val[:-1]
            monthly.append(float(stripped_val))
        else:
            monthly.append(float(val))
    annual_val = data_bits[-1]
    annual = None
    if is_flagged(annual_val):
        stripped_val = annual_val[:-1]
        annual = float(stripped_val)
    else:
        annual = float(annual_val)
    
    return (id, variable, year, monthly, annual)
    
def read_station_string(station_str):
    """Parse a line in a USHCN station list metadata file and convert it
    into a Station object.
    
    Information on these station list files can be found at 
    http://cdiac.ornl.gov/ftp/ushcn_v2_monthly/readme.txt. USHCN rigidly
    formats them, so this method anticipates the USHCN format and reads
    off the data accordingly.
    
    :Param station_str:   
        The USHCN formatted string containing station metadata.
        
    :Return:
        A Station object containing the information from station_str.
         
    """
    station_str = station_str.strip()
    
    coop_id = station_str[0:6]
    lat = float(station_str[7:15])
    lon = float(station_str[16:25])
    elev = float(station_str[26:32])
    state = station_str[33:35]
    name = station_str[36:67].strip()
    
    coop_1 = station_str[67:73]
    coop_2 = station_str[74:80]
    coop_3 = station_str[81:87]
    
    utc_offset = int(station_str[87:90])
    
    new_station = Station(coop_id=coop_id, lon=lon, lat=lat, elev=elev,
                          state=state, name=name, coop_1=coop_1,
                          coop_2=coop_2, coop_3=coop_3,
                          utc_offset=utc_offset)
    return new_station

def write_station_string(station):
    """Write a string containing a station's metadata, as formatted by USHCN.
    
    :Param station:
        The Station object which should be converted to a USHCN-formatted 
        string.
        
    :Return:
        A string containing a station's metadata in the USHCN specified format.
    """    
    station_str_builder = '{coop_id:6} {lat:>8.4f} {lon: >9.4f} {elev:>6.1f} ' \
                          '{state:2} {name:30} {coop_1:6} {coop_2:6} ' \
                          '{coop_3:6} {utc_offset:+d}\n'
    return station_str_builder.format(**station.values_dict)


class ArgumentError(Exception):
    """Exception raised for keyword arguments which are undefined or invalid.
    
    :Ivar key:
        The keyword or variable argument which contained an invalid value.
    :Ivar val:
        The actual invalid value contained in 'key'.
    :Ivar proper:
        (optional) A list of acceptable values for 'key'.
    
    """
    def __init__(self, key, val, proper=None):
        self.key = key
        self.val = val
        self.proper = proper
        
    def __str__(self):
        error_str = ""
        if isinstance(self.val, str):
            error_str = "%s is an invalid value for %s." % (self.val, self.key)
        elif isinstance(self.val, int):
            error_str = "%d is an invalid value for %s." % (self.val, self.key)
        elif isinstance(self.val, float):
            error_str = "%f is an invalid value for %s." % (self.val, self.key)
        else:
            return "An invalid value was provided for %s." % self.key
        
        if self.proper:
            proper_vals_str = ("Accepted values are [" + 
                               ", ".join(self.proper) + "].")
            return error_str + " " + proper_vals_str
        else:
            return error_str
        
    
