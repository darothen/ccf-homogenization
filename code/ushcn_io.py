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
import copy

from ushcn_data import Station, Series

#: The root directory for the USHCN monthly data ftp site
USHCN_FTP_ROOT = "http://cdiac.ornl.gov/ftp/ushcn_v2_monthly/"

#: The name of the USHCN station list file
USHCN_STATION_LIST = "ushcn-stations.txt"

#: A pattern to use for matching and downloading data from the USHCN_FTP_ROOT
#: site.
USHCN_RAW_DATA_PATTERN = "9641C_200912_%s.%s"

#: The element codes identifying what variable is contained in a USHCN dataset,
#: and the scaling factor to convert the recorded values to observations.
ELEM_CODES = { 1: { 'type': 'max', 'scale': 0.1 },
               2: { 'type': 'min', 'scale': 0.1 },
               3: { 'type': 'avg', 'scale': 0.1 },
               4: { 'type': 'pcp', 'scale': 0.01 } }

#: Possible types of USHCN data.
ELEM_TYPES = [data['type'] for data in ELEM_CODES.values()]

#: Possible sources of USHCN data.
SOURCES = ['raw', 'tob', 'F52']

#: Possible flags on USHCN data values.
USHCN_FLAGS = ['E', 'I', 'Q', 'X']

def is_flagged(value, flags=USHCN_FLAGS):
    """Determines if a raw value carries a USHCN flag, signifying that it was
    estimated in some way."""
    for flag in flags:
        if value.endswith(flag):
            return True
    return False

#def get_ushcn_data(source, variable, stations=None):
def get_ushcn_data(params):
    """Download USHCN datasets for use in a homogenization algorithm.
    
    :Param params:
        A parameter object which contains the user-defined parameters for a
        give homogenization project. At a *minimum*, params must contain the
        following valus:
        
        params.data_src
        params.variable
        
        These values define both the source USHCN dataset and the variable to 
        download from the NCDC FTP server.
        
        The source dataset provided by USHCN which the user wishes to download;
        must be either "F52" (bias-adjusted mean monthly values with
        estimates for missing values), "tob" (mean monthly values adjusted only
        for the time of observation bias), or "raw" (unadjusted mean monthly
        values).
        
        The variable whose mean monthly values the user wishes to download; must
        be either "max" (monthly maximum temperatures), "min" (monthly minimum
        temperatures), "avg" (average of monthly maximum and minimum
        temperatures), or "pcp" (total monthly precipitation).
    
        Optionally, a list of station coop_ids can be supplied as params.stations
        if the user only wishes to extract a subset of stations' data.
        
    :Return:
        A collection of USHCN data series.
        
    """
    
    
    if params.data_src not in SOURCES:
        raise ArgumentError("source", params.data_src, proper=SOURCES)
    if params.variable not in ELEM_TYPES:
        raise ArgumentError("variable", params.variable, proper=ELEM_CODES)    
    
    ## SYNTHETIC BENCHMARK BRANCH CONSTANTS INJECTION
    if params.benchmark:
        USHCN_RAW_DATA_PATTERN = "BENCHMARK_%s.%s" 
        USHCN_STATION_LIST = "benchmark-stations.txt"
    
    data_fn = (USHCN_RAW_DATA_PATTERN % (params.data_src, params.variable)) + ".gz"
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
    
    just_stations = copy.deepcopy(all_stations)
        
    ## Now that we have all the stations and their corresponding metadata,
    ## we can process the data series associated with the station.    
    if hasattr(params, 'stations'):
        stations = sorted(params.stations)
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
            station.update_values(series=(station.series+[yearly_data]),
                                  years=(station.years+[year]))
        else:
            if station_id in stations:
                station.update_values(series=[yearly_data], first_year=year,
                                      variable=variable, years=[year])
                appended_stations.append(station_id)
                
    all_series = dict()
    years = range(params.begyr, params.endyr)
    for station_id in stations:
        station = all_stations[station_id]
        station.update_values(subset_years=years)
        series = Series(**station.values_dict)
        all_series[station_id] = series

        #write_series(series, data_dir)
        
    return all_series, just_stations

def read_ushcn_dataset_string(data_str):
    """Parse a line in a USHCN master dataset file and return the information
    about monthly/annual values and station identifying details contained
    in that line for use elsewhere.
    
    USHCN master dataset files are ASCII files which contain an entire dataset
    for a given variable and analysis product on multiple lines. A sample line
    from one of these files is:
    
    01108431930   520    585    577    673    761    805    854    823  \\
        792    644E   584    483E   675E
    
    The first 11 characters, the "identifier", indicate the station where this
    line of data was recorded, the type of data it represents, and the year the
    data was observed. Here, we are looking at station 011084. The 7th
    character, 4, indicates that this is average monthly temperatures. The final
    4 characters, 1930, give the year this entry corresponds to.
    
    The next 12 values are monthly records from January -> December. The final
    13th value is an annual average. If a data entry is "missing," a value of
    -9999 will be recorded in that column. In some cases, values may carry a 
    one-character flag following them.
    
    In the above example, this method will return a list containg the values
    found for (id, variable, year, monthly, annual):
    
    ["011084", 3, 1930, [52.0, 58.5, 57.7, 67.43, 76.1, 80.5, 85.4, 82.3, \\
        79.2, 64.4, 58.4, 48.3], 67.5]
        
    :Param data_str:   
        The string read from a USHCN master dataset file to be parsed.
        
    :Return:
        A list of the values found for (id, variable, year, monthly, annual).
        id will be a 6-digit string, variable and year will be integers, monthly
        will be a 12-element list of floats, and annual will be a float.
        
    TODO: Flag to keep/discard values according to associated flags.
    """
    data_str = data_str.strip()
    data_bits = data_str.split()

    identifier = data_bits[0]
    id = identifier[:6]
    variable = int(identifier[6])
    year = int(identifier[7:11])
    
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
    
    Example:
    
    013160  32.8347  -88.1342   38.1 AL GAINESVILLE LOCK \\
              011694 ------ ------ +6
    (station name allotment split to new line)
    
    This is the station metadata for a station at Gainesville Lock in 
    Alabama with Coop ID 013160. Its (lat, lon) are (32.8347N, 88.1342W) and
    it sits at an eleveation of 38.1 meters above sea level. At some point,
    it was consolidated with a station with Coop ID 011694. Finally, it is
    in the UTC+6 timezone. Passing this string to this method will return
    a Station object with the following data:
    
        Station.coop_id = "013160"
        Station.lat = 32.8347
        Station.lon = -88.1342
        Station.elev = 38.1
        Station.state = "AL"
        Station.name = "GAINESVILLE LOCK"
        Station.coop_1 = "011694"
        Station.coop_2 = Station.coop_3 = "------"
        Station.utc_offset = 6
    
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

def format_station_string(station):
    """Write a string containing a station's metadata, as formatted by USHCN.
    
    :Param station:
        The Station object which should be converted to a USHCN-formatted 
        string.
        
    :Return:
        A string containing a station's metadata in the USHCN specified format.
    """    
    #station_str_builder = '{coop_id:6} {lat:>8.4f} {lon: >9.4f} {elev:>6.1f} ' \
    #                      '{state:2} {name:30} {coop_1:6} {coop_2:6} ' \
    #                      '{coop_3:6} {utc_offset:+d}\n'
    #return station_str_builder.format(**station.values_dict)
    station_str_builder = '{coop_id:6} {lat:>8.4f} {lon: >9.4f} {coop_1:6} {coop_2:6} ' \
                          '{coop_3:6} \n'
    return station_str_builder.format(**station.values_dict)

def write_series(series, out_dir=None):
    """Write a single file similar to a USHCN master dataset file, but which
    contains only the data for the provided series."""
    if not out_dir:
        out_dir = os.getcwd()
    
    coop_id = series.coop_id
    variable_str = series.variable
    out_file_name = "%s_%s.raw" % (coop_id, variable_str)
    out_file = open(os.path.join(out_dir, out_file_name), 'wb')
    
    for (year, data) in zip(series.years, series.series):
        data = map(int, data[:-1])
        data_str = "".join(["{: >5d}  ".format(val) for val in data])
        out_str = "%6s %4d %s\n" % (coop_id, year, data_str)
        out_file.write(out_str)
    out_file.close()    

class ArgumentError(Exception):
    """Exception raised for keyword arguments which are undefined or invalid.
    
    :Ivar key:
        The keyword or variable argument which contained an invalid value.
    :Ivar val:
        The actual invalid value contained in 'key'.
    :Ivar proper:
        (optional) A list of acceptable values for 'key'.
    
    TODO: error_str = "%r is an invalid value for %s" % (self.val, self.key)
    """
    def __init__(self, key, val, proper=None):
        self.key = key
        self.val = val
        self.proper = proper
        
    def __repr__(self):
        error_str = "%r is an invalid value for %s." % (self.val, self.key)
        
        if self.proper:
            proper_vals_str = ("Accepted values are [" + 
                               ", ".join(self.proper) + "].")
            return error_str + " " + proper_vals_str
        else:
            return error_str
        
    
