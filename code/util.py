#!/usr/bin/env python
#
# Daniel Rothenberg, 2011-06-13

"""Utility functions for processing USHCN data.

"""
__docformat__ = "restructuredtext"

# http://docs.python.org/library/math.html
from math import cos, acos, sin, radians, sqrt
# http://docs.python.org/library/operator.html
from operator import itemgetter
# ccf-homogenization imports
from parameters import RADIUS_EARTH

def get_valid_data(data, missing_val=-9999):
    """Return only the valid values in a dataset.
    
    :Param data:
        A list of data values.
    :Param missing_val:
        The placeholder for missing values in the dataset.
    :Return:
        A list of data, with elements matching missing_val removed.
    
    """
    valid_data = [val for val in data if val != missing_val]
    return valid_data

def compute_std(data, missing_val=-9999, valid=False):
    """Computes the unbiased sample standard deviation of a set of data.
    
    :Param data:
        A list of data values.
    :Param missing_val:
        The placeholder for missing values in the dataset.
    :Param valid:
        (optional) Boolean flag indicating that the data has already been 
        sanitized of missing values
    :Return:
        The standard deviation of the dataset. If the dataset has less than 2
        valid entries in it, then return missing_val as the standard deviation
        (we can't compute the standard deviation for 0 or 1 data).
    
    """
    if not valid:
        data = get_valid_data(data, missing_val)
    
    data_mean = compute_mean(data, missing_val, valid=True)
    nval = len(data)
    
    if nval < 2:
        return missing_val    
    
    std = sqrt(sum((d-data_mean)**2 for d in data)/(nval-1))
    return std


def compute_corr(x, y, missing_val=-9999, valid=False):
    """Computes the Pearson Correlation Coefficient between two sets of data.
    
    The Pearson Correlation Coefficient is a measure of the linear dependence
    between two sets of data mapped between [-1, 1]. This method only considers
    values of i where both X[i] and Y[i] are good - that is, not missing.
    
    The code presented here is based in part on a routine written by David
    Jones, http://code.google.com/p/amberfrog/source/browse/trunk/zontem/code/correlation.py.
    
    :Param x,y:
        The datasets for which the correlation should be computed, supplied as a
        list of floats or ints.
    :Param missing_val:
        The placeholder for missing values in either dataset.
    :Param valid:
        (optional) Boolean flag indicating that both datasets have already been 
        sanitized of missing values
    :Return:
        Correlation coefficient (r in [-1.0, 1.0]) if both x and y have 
        equal amounts of valid data (more than 0); otherwise, returns None.
    
    """
    ## Align x and y so that we are computing correlations only where
    ## both x[i] and y[i] are valid data. This will also eliminate all the
    ## missing values in the end x and y lists, so we can pass the "valid"
    ## argument on to later computations.
    if not valid:
        aligned_zip = [(xi, yi) for xi, yi in zip(x, y) if ((xi != missing_val) and 
                                                        (yi != missing_val))]
        x, y = zip(*aligned_zip)
    
    n = len(x)
    assert len(y) == n # Computation assumes len(x_valid) == len(y_valid)
    
    ## If there were fewer than 2 valid data values in each set, then we can't
    ## compute the standard deviation and therefore can't compute the
    ## correlation coefficient. The USHCN code has a higher threshold, and
    ## further requires that there actually be 8 months of overlap for stations
    ## to be correlated as pairs. We'll use that higher threshold here.
    if n < 2:
        return None
        
    ## Now, there are no missing data in x_valid or y_valid, so we can pass
    ## a flag to the mean and std functions to avoid having to filter through
    ## the data second and third times.
    x_bar = compute_mean(x, missing_val, valid=True)
    y_bar = compute_mean(y, missing_val, valid=True)
    
    numerator = sum((xi-x_bar)*(yi-y_bar) for (xi, yi) in zip(x, y))    
    
    x_std = compute_std(x, missing_val, valid=True)
    y_std = compute_std(y, missing_val, valid=True)
    
    rank = numerator/((n-1)*x_std*y_std) # Divide-by-zero is possible, and
                                         # should throw an exception.
    return rank
    
def compute_mean(data, missing_val=-9999, valid=False):
    """Computes the mean of a given set of data, with the possibility that
    the dataset contains missing values.
    
    :Param data:
        The data over which to compute the mean.
    :Param missing_val:
        The placeholder for missing values in the dataset.
    :Param valid:
        (optional) Boolean flag indicating that the data has already been 
        sanitized of missing values
    :Return:
        The mean of the dataset and the number of values used to compute it. If
        all the data was missing, will return missing_val.
    
    """
    if not valid:
        data = get_valid_data(data, missing_val)
    
    total = sum(data)*1.0
    nval = len(data)*1.0 # Convert to float just to avoid integer truncation
    
    ## If there aren't any valid_data, then we'll return missing_val now
    if not nval: 
        return missing_val
    
    mean = total/nval
    return mean 
        
def compute_first_diff(monthly_data, missing_val=-9999):
    """Computes the first-order timeseries differences for a dataset.
    
    Given a dataset comprised as a list of datavalues which possibly contains
    a placeholder "missing" value specified by the user, computes the first
    order difference timeseries from that dataset. This timeseries is defined
    as:
    
    F[t] = X[t] - X[t-1]
    
    where t runs from 1 through the length of X. If either X[t] or X[t-1] is
    a missing value, then F[t] is defined to be missing_val.
    
    :Param monthly_data:
        The monthly data to compute anomalies for.
    :param missing_val:
        The placeholder for missing data.
    :Return:
        A list of dimensions (len(monthly_data)-1) containing the first-order
        difference timeseries from monthly_data.
    
    """
    
    data_left = monthly_data[0:-1]
    data_right = monthly_data[1:]
    
    first_diffs = []
    for (left, right) in zip(data_left, data_right):
        if (left != missing_val and right != missing_val):
            ## NORMAL, TEXTBOOK FIRST DIFFERENCE FORMULA
            #first_diffs.append(right-left)
            ## FILTERED FIRST DIFFERENCE FORMULA USED IN 
            ## ushcn_corr_2004.v3.f, subroutine frstdif
            first_diffs.append((right-left)/2.)
        else:
            first_diffs.append(missing_val)
            
    return first_diffs
    

def compute_monthly_anomalies(yearly_data, missing_val=-9999):
    """Computes monthly anomalies given a series of data.
    
    Monthly anomalies can be ambiguous defined, but here, a monthly anomaly
    is computed by taking the mean value for all of a specific month in a
    dataset. That is, if you have a list of yearly/monthly data, the anomaly for 
    January will be computed by taking the mean value of all January data
    entries, and subtracting each January value from that mean.
    
    :Param yearly_data:
        The yearly data to compute anomalies for, with dimensions (nyears x 12)
    :param missing_val:
        The placeholder for missing data which shouldn't be accumulated.
    :Return:
        A list of the same dimensions as yearly_data, but containing
        anomalies instead of the original values.
    
    -----------------
    
    NOTE - This is probably terribly inefficient and kind of silly, but works
    exactly as it is supposed to. I'll come back and try to write a clearer
    function in the future (this is where numpy and being able to slice multi-
    dimensional arrays would be super handy!)
    
    """
    
    nyears = len(yearly_data)
    nmonths = 12
    anomalies = []
    for i_month in range(nmonths):
        all_months_data = map(itemgetter(i_month), yearly_data)
        month_mean = compute_mean(all_months_data, missing_val)
        if i_month == 0:
            for i_year in range(nyears):
                anomalies.append([anomaly(all_months_data[i_year],
                                          month_mean, missing_val)])
        else:
            for i_year in range(nyears):
                anomalies[i_year].append(anomaly(all_months_data[i_year],
                                                  month_mean, missing_val))
    return anomalies

def anomaly(val, base, missing_val=-9999):
    """Computes an anomaly given a value and a baseline.
    
    :Param val:
        The value to compute the anomaly from.
    :Param base:
        The base used to compute the anomaly.
    :Return:
        The anomaly, val-base. If val==missing_val, will return missing_val
        instead.
    
    """
    if val == missing_val:
        return missing_val
    else:
        return (val - base)

def compute_arc_dist(station1=None, station2=None,
                     lat1=None, lon1=None, lat2=None, lon2=None):
    """Compute the arc-distance between two stations.
    
    Given either two Stations or their latitude and longitude in
    degrees, computes the distance along the sphere between the
    two Stations. This function assumes that latitudes are negative
    in the Southern hemisphere and positive in the North, and that
    longitudes are negative in the Western hemisphere and positive
    in the East.
    
    If both a set of two Stations or set of four coordinates is
    supplied, the function will use the four coordinates in lieu
    of the Stations' metadata.
    
    The formula used here derives from computing the angle between
    two vectors (one extending from the center of the earth to each
    station, assuming a spherical Earth) by using the definition
    of the dot product ( dot(p, q) = mag(p)*mag(q)*cos(theta_p,q) ),
    and using that angle to compute an arc distance.
        
    :Param station1, station2:
        The Station objects which contain metadata about observation
        sites, including latitude and longitude in degrees.
    :Param lat1, lon1:
        The latitude and longitude of the first station, in degrees.
    :Param lat2, lon2: 
        The latitude and longitude of the second staiton, in degrees.
    :Return:
        The arc-distance between the two Stations, in kilometers.
    
    """
    if not (lat1 and lat2 and lon1 and lon2):
        lat1, lon1 = (station1.lat, station1.lon)
        lat2, lon2 = (station2.lat, station2.lon)
        
    lat1rad, lon1rad, lat2rad, lon2rad = map(radians,
                                             [lat1, lon1, lat2, lon2])
    
    x_diff = cos(lon1rad)*cos(lat1rad)*cos(lon2rad)*cos(lat2rad)
    y_diff = sin(lon1rad)*cos(lat1rad)*sin(lon2rad)*cos(lat2rad)
    z_diff = sin(lat1rad)*sin(lat2rad)
    
    arc_angle = acos(x_diff + y_diff + z_diff)
    arc_dist = arc_angle*RADIUS_EARTH
    
    return arc_dist
    
    
