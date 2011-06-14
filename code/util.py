#!/usr/bin/env python
#
# Daniel Rothenberg, 2011-06-13

"""Utility functions for processing USHCN data.

"""
__docformat__ = "restructuredtext"

# http://docs.python.org/library/math.html
from math import cos, acos, sin, radians, sqrt
# ccf-homogenization imports
from parameters import RADIUS_EARTH

def compute_corr(x, y, missing_val=-9999):
    """Computes the Pearson Correlation Coefficient between two sets of data.
    
    The Pearson Correlation Coefficient is a measure of the linear dependence
    between two sets of data mapped between [-1, 1]. This method only considers
    values of i where both X[i] and Y[i] are good - that is, not missing.
    
    :Param x,y:
        The datasets for which the correlation should be computed.
    :Param missing_val:
        The placeholder for missing values in either dataset.
    :Return:
        Correlation coefficient, standard deviation of x, standard deviation of
        y.
    
    """
    ## Need a very small value because the numerics here can be contaminated
    ## by round-off error, which can cause trouble if we're dividing by small
    ## numbers.
    eps = 1e-6
    
    x_good = []
    y_good = []
    for (x_val, y_val) in zip(x,y):
        if (x_val != missing_val and y_val != missing_val):
            x_good.append(x_val)
            y_good.append(y_val)
    
    ## Is there actually data to work with? If not, don't bother doing
    ## any math, and return r = 0, x_std = 0, y_std = 0
    if not (x_good and y_good):
        return 0, 0, 0
    
    x_bar, nx = compute_mean(x_good, missing_val)
    y_bar, ny = compute_mean(y_good, missing_val)
    
    ## Compute standard deviation of x, y
    x_std = 0.0
    y_std = 0.0
    n_std = 0.
    for (x_val, y_val) in zip(x_good, y_good):
        x_std = x_std + (x_val-x_bar)**2
        y_std = y_std + (y_val-y_bar)**2
        n_std = n_std+1
    x_std = sqrt(x_std/n_std)
    y_std = sqrt(y_std/n_std)
    
    ## Use the standard deviations and means to now compute the actual
    ## correlation coefficient.
    sum = 0.0
    n_sum = 0.
    for (x_val, y_val) in zip(x_good, y_good):
        sum = sum + ((x_val-x_bar)*(y_val-y_bar))
        n_sum = n_sum+1
    ## Replace standard deviations if they're small enough to cause numerical
    ## troubles.
    if x_std < eps:
        x_std = eps
    if y_std < eps:
        y_std = eps    
    r = sum/((n_sum-1)*x_std*y_std)
    
    return r, x_std, y_std

def compute_mean(data, missing_val=-9999):
    """Computes the mean of a given set of data, with the possibility that
    the dataset contains missing values.
    
    :Param data:
        The data over which to compute the mean.
    :Param missing_val:
        The placeholder for missing values in the dataset.
    :Return:
        The mean of the dataset and the number of values used to compute it. If
        all the data was missing, will return missing_val.
    
    """
    
    total = 0.
    nval = 0.
    for val in data:
        if val != missing_val:
            total = total + val
            nval = nval + 1.
    
    if nval > 0:
        mean = total/nval
    else:
        mean = missing_val
        
    return (mean, nval)
        

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
    

def compute_monthly_anomalies(monthly_data, missing_val=-9999):
    """Computes monthly anomalies given a monthly series of data.
    
    :Param monthly_data:
        The monthly data to compute anomalies for.
    :param missing_val:
        The placeholder for missing data which shouldn't be accumulated.
    :Return:
        A list of the same dimensions as monthly_data, but containing
    
    """
    mean, nval = compute_mean(monthly_data, missing_val)
            
    ## Did we actually accumulate values? If not, then all the data
    ## is missing values, so just return the original list.
    if not nval:
        return monthly_data
    
    anomalies = []
    for val in monthly_data:
        if val != missing_val:
            anomalies.append(val - mean)
        else:
            anomalies.append(missing_val)
    
    return anomalies    
    

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
    
    
