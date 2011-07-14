#!/usr/bin/env python
#
# Copyright (C) 2011 Daniel Rothenberg.
# See Google Code project page for license, 
# http://code.google.com/p/ccf-homogenization/

"""Port of functions necessary for identifying and classifying changepoint
model types during splitmerge process

"""
__docformat__ = "restructuredtext"

# http://docs.python.org/library/math.html
from math import log10, sqrt

# ccf-homogenization imports
from util import compute_mean, get_valid_data, compute_std

class Stats(dict):
    """Base object for returning results of statistical tests
    
    """
    def __init__(self, test_name, **params):
        self.test_name = test_name
        self.__dict__.update(params)
        
    def __repr__(self):
        return "%s (%r)" % (self.test_name, self.__dict__)
    

def least_squares(x, y, missing_val=-9999):
    """Least-squares regression using Ordinary Least Squares for two 1D data
    sets.
    
    Given two sets of scalar data, performs a simple least squares
    regression over all pairs of data where neither x_i nor y_i are equal
    to the supplied missing value. Assumes that x and y are the same length.
    
    :Params x, y:
        The predictor and predictand data, respectively. Each should be a one
        dimensional list 
    :Param missing_val:
        The placeholder for missing values in the dataset.
    :Return:
        A 'Stats' object with the attributes x_mean, y_mean, slope, y_int,
        sseflat, sseslope, and nval.
    
    """        
    assert len(x) == len(y)
    
    # fcn to determine whether a value in a pair is invalid
    valid_pair = lambda pair: ((pair[0] != missing_val) and 
                               (pair[1] != missing_val))
    
    pairs = zip(x, y)
    # Filter out data pairs where either x_i or y_i equal missing_val
    good_data = [pair for pair in pairs if valid_pair(pair)]
    good_x, good_y = zip(*good_data)
    
    # Compute means in both x, y
    nval = len(good_data)
    x_mean = compute_mean(good_x, valid=True)
    y_mean = compute_mean(good_y, valid=True)
    
    '''
    nval = 0.0
    sum_x, sum_y = 0.0, 0.0
    for (x_val, y_val) in zip(x, y):
        if x_val != missing_val and y_val != missing_val:
            nval = nval + 1
            sum_x = sum_x + x_val
            sum_y = sum_y + y_val
            
    x_mean = sum_x/nval
    y_mean = sum_y/nval
    '''
                
    # calculate slope and y intercept   
    numer = 0.0
    denom = 0.0
    for (x_val, y_val) in good_data:
        numer = numer + (x_val-x_mean)*(y_val-y_mean)
        denom = denom + (x_val-x_mean)*(x_val-x_mean)
    '''
    for (x_val, y_val) in zip(x, y):
        if x_val != missing_val and y_val != missing_val:
            rnum1 = rnum1 + (x_val-x_mean)*(y_val-y_mean)
            rden1 = rden1 + (x_val-x_mean)*(x_val-x_mean)
    '''
    
    slope = numer/denom
    y_int = y_mean - slope*x_mean
    
    # calculate residuals
    sseflat = 0.0
    sseslope = 0.0
    for (x_val, y_val) in good_data:
        sseflat = sseflat + (y_mean-y_val)**2
        resid = (y_int + slope*x_val) - y_val
        sseslope = sseslope + resid**2
    '''
    for (x_val, y_val) in zip(x, y):
        if x_val != missing_val and y_val != missing_val:
            sseflat = sseflat + (y_mean - y_val)**2
            resid = (y_int + slope*x_val) - y_val
            sseslope = sseslope + resid**2
    '''
    
    # return a Stats object with all of these computed values
    return Stats(test_name="least-squares", x_mean=x_mean, y_mean=y_mean,
                 slope=slope, y_int=y_int, sseflat=sseflat, sseslope=sseslope,
                  nval=nval)
    
def kth_line(x, y, missing_val=-9999):
    """Estimates a linear regression using the Theil-Sen estimator (otherwise
    known as the Kendall-Theill Robust Line fit).
    
    Given two sets of data, computes the Kendall-Theill Robut Line fit for all
    pairs of data where both values are non-missing. Assumes that both given
    data sets are the same length.
    
    :Params x, y:
        The predictor and predictand data, respectively. Each should be a one
        dimensional list.
    :Param missing_val:
        The placeholder for missing values in the dataset.
    :Return:
        A 'Stats' object with the attributes x_med, y_med, slope, y_int,
        sseflat, sseslope, and nval.
    
    """
    assert len(x) == len(y)
    
    # fcn to determine whether a value in a pair is invalid
    valid_pair = lambda pair: ((pair[0] != missing_val) and 
                               (pair[1] != missing_val))
    
    pairs = zip(x, y)
    # Filter out data pairs where either x_i or y_i equal missing_val
    good_data = [pair for pair in pairs if valid_pair(pair)]
    good_x, good_y = zip(*good_data)
    nval = len(good_data)
    
    '''
    rx, ry = [], []
    for (x_val, y_val) in zip(x, y):
        if x_val != missing_val and y_val != missing_val:
            rx.append(x_val)
            ry.append(y_val)
    nval = len(rx)
    '''
    
    ## calculate slope
    # 1) if less than 60 values, use sort-ascending technique
    #if nval < 60:
    if nval > 0:
        nslp = 0
        r_temp = []
        for i in range(nval-1):
            for j in range(i, nval):
                #if rx[i] != rx[j]:
                if good_x[i] != good_x[j]:
                    nslp = nslp + 1
                    #r_temp.append( (ry[j]-ry[i])/(rx[j]-rx[i]) )
                    r_temp.append( (good_y[j]-good_y[i])/(good_x[j]-good_x[i]) )
        # sort into ascending order
        r_temp = sorted(r_temp)
        # find the median value -
        if (nslp%2) == 1:
            imed = nslp/2
            slope = r_temp[imed]
        else:
            imed = (nslp-1)/2
            slope = (r_temp[imed]+r_temp[imed+1])/2.0
    # 2) if more than 60 values, use binned distribution
    #else:


    # calculate y-intercept; use original lists so we can get true median
    # for all of the input data
    rx = sorted(good_x)
    ry = sorted(good_y)
    if (nval%2) == 1:
        imed = nval/2
        x_med = rx[imed]
        y_med = ry[imed]
    else:
        imed = (nval-1)/2
        x_med = (rx[imed]+rx[imed+1])/2.0
        y_med = (ry[imed]+ry[imed+1])/2.0
    
    y_int = y_med - slope*x_med
    
    # Calculate residuals
    sseflat = 0.0
    sseslope = 0.0
    #for (x_val, y_val) in zip(x, y):
    #    if x_val != missing_val and y_val != missing_val:
    for (x_val, y_val) in good_data:
        sseflat = sseflat + (y_med - y_val)**2
        resid = (y_int + slope*x_val) - y_val
        sseslope = sseslope + resid**2
    
    # return a Stats object with these computed values
    return Stats(test_name="kth-line", x_med=x_med, y_med=y_med, slope=slope,
                 y_int=y_int, sseflat=sseflat, sseslope=sseslope, nval=nval)  
      
def bayes(nval, sum_square_resids, dof):
    """Computes the Bayes Information Criterion.
    
    """    
    nval = float(nval)
    rsq1 = nval*log10(sum_square_resids/nval)
    rsq2 = dof*log10(nval)
    
    return (rsq1+rsq2, rsq1, rsq2)

def t_test(x, y, missing_val=-9999):
    
    valid_x = get_valid_data(x, missing_val)
    valid_y = get_valid_data(y, missing_val)
    
    # Compute sample size, standard deviations, and means for
    # each of the two populations
    mean_x = compute_mean(valid_x, valid=True)
    mean_y = compute_mean(valid_y, valid=True)
    
    std_x = compute_std(valid_x, valid=True)
    std_y = compute_std(valid_y, valid=True)
    
    num_x = len(valid_x)
    num_y = len(valid_y)
    
    # Now compute the t-statistic
    numer = abs(mean_x - mean_y)
    denom = sqrt(((std_x**2)/num_x)+((std_y**2)/num_y))
    
    t_val = numer/denom
    
    # return a Stats object with these computed values
    return Stats(test_name="t-test", mean_x=mean_x, mean_y=mean_y, std_x=std_x,
                 std_y=std_y, num_x=num_x, num_y=num_y, t_val=t_val)
    
    
def lookup_critical(n, stat="t"):
    """Uses lookup tables and linear interpolation to compute the critical value
    at 95% significance for a given test statistic and number of observations.
    
    :Param n:
        The basis of the test statistic critical value. In most cases, this
        should equal the number of test realizations minus the degrees of
        freedom or number of parameters necessary to fit the test sample
        population.
    :Param stat:
        The name of the test statistic to use. Possible values are -
            "t" - 2-tailed t_test
            "f1" - 
    :Return:
        The critical test_statistic value.            
    
    """
    lookup_tables = {
                     "t": [12.71, 4.30, 3.18, 2.78, 2.57, 2.45, 2.37,2.31,
                           2.26, 2.23, 2.20, 2.18, 2.16, 2.15, 2.13, 2.12, 
                           2.11, 2.10, 2.09, 2.09, 2.08, 2.07, 2.07, 2.06, 
                           2.06, 2.06, 2.05, 2.05, 2.05, 2.04, 2.02, 2.00, 
                           1.98, 1.96], 
                     "f1": [161.4, 18.50, 10.1, 7.7, 6.6, 6.0, 5.6, 5.3,
                              5.1,   5.0,  4.8, 4.7, 4.7, 4.6, 4.5, 4.5,
                              4.5,   4.4,  4.4, 4.4, 4.3, 4.3, 4.3, 4.3,
                              4.2,   4.2,  4.2, 4.2, 4.2, 4.2, 4.1, 4.0,
                              3.9,   3.8],
                     "f2": [199.5, 19.00, 9.55, 6.94, 5.79, 5.14, 4.74, 4.46, 
                             4.26,  4.10, 3.98, 3.89, 3.81, 3.74, 3.68, 3.63,
                             3.59,  3.55, 3.52, 3.49, 3.47, 3.44, 3.42, 3.40,
                             3.39,  3.37, 3.35, 3.34, 3.33, 3.32, 3.23, 3.15,  
                             3.07, 3.00],                     
                     }
    thresholds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                  19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 40, 60, 120,
                  9999]
    crit_vals = lookup_tables[stat]
    
    # Can we short-circuit this lookup? Yes we can, under two conditions -
    #     1) n > 9999; in this case, we use the final critical value in the
    #        lookup table.
    #     2) n <= 30; in this case, we already know the critical
    #        value and don't need to interpolate.
    if n > 9999:
        return crit_vals[-1]
    elif n <= 30:
        return crit_vals[n-1]
    else:
    # Else, we need to do a linear interpolation.
    # Use the thresholds list to find the index to use in the lookup table
        for i in range(len(thresholds)):
            if n < thresholds[i]:
                break
        
        thresh_left, thresh_right = thresholds[i-1], thresholds[i]
        crit_left, crit_right = crit_vals[i-1], crit_vals[i]
        slope = ((crit_right-crit_left)/(thresh_right-thresh_left))
        
        return (crit_left + slope*(n-thresh_left))

    

    
        
    