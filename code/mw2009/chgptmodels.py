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
# http://docs.python.org/library/operator.html
from operator import itemgetter

# ccf-homogenization imports
from util import compute_mean, get_valid_data, compute_std

# Possible models to fit to changepoints
minbic_models = [
                 # Homogeneous linear types -
                 'KTHSLR0', # simple flat line
                 'KTHSLR1', # sloped straight line
                 # Two phase regression types - 
                 'KTHTPR0', # amplitude only, 0 sloped segments (step change)
                 'KTHTPR1', # amplitude shift with equal sloped segments
                 'KTHTPR2', # amplitude shift with non-equal sloped segments
                 'KTHTPR3', # amplitude shift with flat-to-slope segments
                 'KTHTPR4', # amplitude shift with slope-to-flat segments  
                 ]

class Stats(dict):
    """Base object for returning results of statistical tests
    
    """
    def __init__(self, test_name, **params):
        dict.__init__(self, **params)
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
    
    Given two sets of data, computes the Kendall-Theill Robust Line fit for all
    pairs of data where both values are non-missing. Assumes that both given
    data sets are the same length.
    
    :Params x, y:
        The predictor and predictand data, respectively. Each should be a one
        dimensional list.
    :Param missing_val:
        The placeholder for missing values in the dataset.
    :Return:
        A 'Stats' object with the attributes 
            x_med - median of predictor
            y_med - median of predictand
            y_int, slope - parameters of regression fit to x, y data
            sseflat - sum squared error computed between y_med and all y vals
            sseslope - sum squared error of the residuals fit to this regression
            nval - the number of valid pairs of values found
    
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
        
    # calculate paired slopes, find median value
    if nval > 0:
        nslp = 0
        slopes = []
        for i in range(nval-1):
            for j in range(i, nval):
                if good_x[i] != good_x[j]:
                    nslp = nslp + 1
                    slopes.append( (good_y[j]-good_y[i])/(good_x[j]-good_x[i]) )
        # sort into ascending order
        slopes = sorted(slopes)
        # find the median value -
        if (nslp%2) == 1:
            imed = nslp/2
            slope = slopes[imed]
        else:
            imed = (nslp-1)/2
            slope = (slopes[imed]+slopes[imed+1])/2.0

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
    for (x_val, y_val) in good_data:
        sseflat = sseflat + (y_med - y_val)**2
        resid = (y_int + slope*x_val) - y_val
        sseslope = sseslope + resid**2
    
    # return a Stats object with these computed values
    return Stats(test_name="kth-line", x_med=x_med, y_med=y_med, slope=slope,
                 y_int=y_int, sseflat=sseflat, sseslope=sseslope, nval=nval)  
      
def bayes(nval, sse, dof):
    """Computes the Bayes Information Criterion.
    
    :Param nval:
        Number of data points in the segment being fit.
    :Param sum_square_resids:
        The sum of the squared errors of the fit.
    :Param dof:
        The degrees of freedom in the fit - should be equal to the number of
        parameters being fit.
    :Return:
        A tuple containing the BIC value, the BIC contribution from the error of
        the model fit, and the BIC penalty due to the number of parameters fit
        in the model.
    
    """    
    nval = float(nval)
    rsq1 = nval*log10(sse/nval) # BIC contribution from error
    rsq2 = dof*log10(nval) # BIC penalty for parameters fit
    
    return (rsq1+rsq2, rsq1, rsq2)

def t_test(x, y, missing_val=-9999):
    """Performs a Student's T-Test on a set of data.
    
    :Params x, y:
        The two datasets being compared in the statistical test. Should be lists
        of floats, but some values can equal missing_val.
    :Param missing_val:
        (optional) The placeholder value for missing values which should not be
        considered in this analysis; default is -9999.
    :Return:
        A 'Stats' object containing the results of this test as attributes -
            mean_x, std_x, num_y - the mean, standard deviation, and number of
                valid values found in the x dataset
            mean_y, std_y, num_y - the mean, standard deviation, and number of
                valid values found in the y dataset
            t_val - the computed t value test statistic for this data
    
    """
    
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
            "f1" - right-tailed f-test, 1 degree of freedom in denominator
            "f2" - right-tailed f-test, 2 degrees of freedom in denominator
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
    
    # Test that we received a valid n and stat value
    if n <= 1:
        raise ValueError("'n' must be >= 1")
    if not stat in lookup_tables:
        raise ValueError("'stat' must be one of %r" % lookup_tables.keys())
    
    crit_vals = lookup_tables[stat]
    
    # Can we short-circuit this lookup? Yes we can, under two conditions -
    #     1) n > 9999; in this case, we use the final critical value in the
    #        lookup table, which is the "n -> infinity" case.
    #     2) 0 < int(n) <= 30; in this case, we already know the critical
    #        value (it's in our table) and don't need to interpolate.
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

def minbic(x, y, bp_index, missing_val=-9999):
    """Performs a changepoint-detection analysis to classify a suspected
    changepoint as one of a pre-specified list of models (minbic_models, above).
    
    This method performs several statistical tests in an attempt to fit a
    model to the suspected changepoint in the data passed to it. In this version
    of the method, the statistical tests are based on the Kendall-Theill method,
    although future versions will support ordinary least squares and alternative
    model-fit analyses. 
    
    A key feature of the method is that the tests based on the Kendall-Theill
    statistic are allowed to cascade; computations in the earlier tests are
    re-used (such as sum square error on the residuals for the line fit to the
    entire data segment) ot minimize the amount of computations necessary here.
    These re-used results are tracked in comments throughout the method. Also,
    the tests used to fit certain models are described in comments in this 
    method as they are encountered.
    
    Upon completion of the changepoint tests, this method will identify which
    test minimizes the Bayesian Information Criterion, based on the error
    computed for each test. It will then return that changepoint's statistics
    for later use.
    
    :Params x, y:
        The predictor and predictand dataset where there is a suspected
        changepoint. Must be equal length lists of floats, but may contain
        missing_val in place of missing values.
    :Param bp_index:
        The index of the suspected changepoint in the x and y datasets.
    :Param missing_val:
        (optional) The placeholder for missing values in the datasets.
    :Return:
    
    """    
    
    assert len(y) == len(x)

    ## Print header for BIC changepoint testing -
    left_header = " QTYP     QVAL    QRSE     QPF     MU1     MU2  ALPHA1"
    right_header = "  ALPHA2   MSTAT   MCRIT    MOFF KNT1 KNT2"
    print (left_header+right_header)
    
    # For storing the results of these tests.
    changepoint_dict = dict()
    
    ############################################################################
    ## TEST 1 - KTHSLR1, Sloped straight line
    ##    Method - Use Kendall-Theill method to compute a linear regression to
    ##            all of the data. We tacitly assume this model as the null
    ##            hypothesis; we must reject this model in order to accept one
    ##            of the more-complicated ones to follow.    
    cmodel = "KTHSLR1"
    kth_regression = kth_line(x, y, missing_val)
    
    num_obs =  len(x) # total number of observations
    n_vals = kth_regression.nval # number of non-missing observations
    slope = kth_regression.slope
    y_int = kth_regression.y_int
    sse_resid = kth_regression.sseslope
    amp_est = 0.0
    
    # Bayes Info Criterion - total, due to error in variance, and penalty
    b_tot, b_err, b_pen = bayes(num_obs, sse_resid, 2)
    # output string
    head = "%7s %6.2f %7.2f %7.2f" % (cmodel, b_tot, b_err, b_pen)
    stats =  " %7.2f ------- %7.3f ------- ------- -------" % (y_int, slope)
    tail = " %7.2f %5d ----" % (amp_est, n_vals)
    print (head+stats+tail)
    
    ## Record the results of this test
    changepoint_dict[cmodel] = dict(bic=b_tot, 
                                    mu=[y_int, y_int],
                                    alpha=[slope, slope], 
                                    sse_bic=(sse_resid/n_vals),
                                    seg_lens=[n_vals, 0],
                                    test_stat=0.0,
                                    crit_val=0.0,
                                    amp_est=amp_est) 
    
    ############################################################################  
    ## TEST 2 - KTHTPR0, Amplitude change, no slope change
    ##    Method - Use Kendall-Theill method to compute a linear regression to
    ##            the segments on either side of the suspected break point. Then,
    ##            applies a student's t-test to see whether or not the mean in
    ##            either segment is different from the other.
    cmodel = "KTHTPR0"
    
    # Divide the segment before/after the suspect breakpoint
    left_x = x[:bp_index+1]
    left_y = y[:bp_index+1]
    right_x = x[bp_index+1:]
    right_y = y[bp_index+1:]
    
    # Linear regression on each segment to find median values to estimate
    # step change
    kthl_left = kth_line(left_x, left_y, missing_val)
    kthl_right = kth_line(right_x, right_y, missing_val)
    
    left_y_med = kthl_left.y_med
    right_y_med = kthl_right.y_med
    amp_est = left_y_med - right_y_med
    
    n_left = kthl_left.nval
    n_right = kthl_right.nval
    n_total = n_left + n_right
    
    # T-test on the mean of left, right segments
    stat_test = t_test(left_y, right_y, missing_val)
    t_val = stat_test.t_val
    t_crit = lookup_critical(n_total-2, "t")
    
    # Use error about the median on left, right segments
    sse_sum = kthl_left.sseflat + kthl_right.sseflat
    b_tot, b_err, b_pen = bayes(n_total, sse_sum, 3)
    # output string
    head = "%7s %6.2f %7.2f %7.2f" % (cmodel, b_tot, b_err, b_pen)
    stats =  " %7.2f %7.2f ------- ------- %7.2f %7.2f" % (left_y_med, right_y_med,
                                                          t_val, t_crit)
    tail = "% 7.2f %5d %4d" % (amp_est, n_left, n_right)
    print (head+stats+tail)
    
    ## add to results dictionary
    changepoint_dict[cmodel] = dict(bic=b_tot, 
                                    mu=[left_y_med, right_y_med],
                                    alpha=[0.0, 0.0],
                                    sse_bic=sqrt(sse_sum/n_total),
                                    seg_lens=[n_left, n_right],
                                    test_stat=t_val,
                                    crit_val=t_crit,
                                    amp_est=amp_est)
    
    ############################################################################  
    ## TEST 3 - KTHTPR1, Amplitude change, equal slope segments
    ##    Method - Use Kendall-Theill method to compute the slope for the entire
    ##            data segment, and use that but coupled with independent 
    ##            y-intercepts for the left/right segments to find residuals
    ##            and error terms.
    ##
    ##    Note - This part of the code does a lot of heavy lifting which comes
    ##            in handy when analyzing the subsequent changepoint models. 
    ##            The various computations are documented here; the next models
    ##            will indicate when they borrow results from this segment of
    ##            the code.
    cmodel = "KTHTPR1"
    
    all_data = y
    num_obs = len(y)
    left_data = y[:bp_index+1]
    right_data = y[bp_index+1:]
    
    ## 1) Compute the mean for *all* of the data
    all_valid_data = get_valid_data(all_data, missing_val)
    all_mean = compute_mean(all_valid_data, valid=True)
    
    ## 2) use Kendall-Theill method with single slope
    # This method is slightly different than the kth_line() method above.
    # First, we get only the valid data, and we pair it with the natural ordering
    # of the data, i.e. 1, 2, 3...
    valid_all = all_valid_data
    n_all = len(valid_all)
    range_all = range(1, n_all+1)
    
    valid_left = get_valid_data(left_data, missing_val)
    valid_right = get_valid_data(right_data, missing_val)
    
    n_left, n_right = len(valid_left), len(valid_right)
    range_left = range(1, n_left+1)
    range_right = range(n_left+1, n_all+1)
    
    # Second, generate paired slopes for the first segment.
    nslp = 0
    slopes = []
    for i in range(n_left-1):
        for j in range(i, n_left):
            if range_left[j] != range_left[i]:
                nslp = nslp + 1
                slopes.append( (valid_left[j]-valid_left[i])/
                               (range_left[j]-range_left[i]) )
    # Third, generate paired slopes for the second segment.
    for i in range(n_right-1):
        # BUG: MW2009 code in chgptmodels.kendallthiell, line 2229 starts the 'j'
        #     index at ibeg2+1. This corresponds to 1 here. Above in the first
        #     segment and in kth_line(), it starts the 'j' index right where 'i' 
        #     left off.
        for j in range(1, n_right):
            if range_right[j] != range_right[i]:
                nslp = nslp + 1
                slopes.append( (valid_right[j]-valid_right[i])/
                               (range_right[j]-range_right[i]) )
                
    #Fourth, find the median slope from all the ones we computed
    slope_ind = 1
    if not slope_ind:
        slope_est = 0.0
    else:
        slopes = sorted(slopes)
        
        # BUG: This is the median-index finding method used in the original
        #    PHA. It's different from what is used in kth_line() (which is
        #    the correct way to find the median in a list of data - average the
        #    middle two values if the list length is even).
        imed = (nslp - 1)/2
        if (nslp%2)==1: imed = imed+1 # offset by one to right if odd
        
        slope_est = slopes[imed]
        
    print "slope, ic, imet: %7.2f %5d %5d" % (slope_est, nslp, imed)
    
    # Fifth, compute the first segment intercept, y-median - slope*x-median
    # BUG: Another instance where we do not handle the case of an even number of
    #    values when computing the median
    imed = (n_left - 1)/2
    if (n_left%2)==1: imed = imed+1
    range_med = range_left[imed]
    valid_left = sorted(valid_left)
    data_med = valid_left[imed]
    left_y_int = data_med-slope_est*range_med
    print "Seg1 - Xmed, Ymed, slope, Yint: %7.2f %7.2f %7.2f %7.3f" % (range_med, 
                                                                       data_med,
                                                                       slope_est,
                                                                       left_y_int)
    # BUG: Again in chgptmodel.kendalltheill(), there is a bug on line 2339. Starting
    #     here, we over-write the medians we found in both lists, and use the second
    #     segment for all our computations! I reproduce that behavior here by using the
    #     generic range_med and data_med values for rXmed and rYmed. This needs 
    #     to be fixed by computing residuals in each segment with respect to the
    #     median from that segment.
    # Sixth, compute the second segment intercept
    # BUG: Another instance where we do not handle the case of an even number of
    #    values when computing the median
    imed = (n_right - 2)/2
    if (n_right%2)==1: imed = imed + 1
    
    range_med = range_right[imed]
    valid_right = sorted(valid_right)
    data_med = valid_right[imed]
    right_y_int = data_med-slope_est*range_med
    print "Seg2 - Xmed, Ymed, slope, Yint: %7.2f %7.2f %7.2f %7.3f" % (range_med, 
                                                                       data_med,
                                                                       slope_est,
                                                                       right_y_int)
    
    # Seventh, we compute root mean square error of the residuals
    residuals = [missing_val]*n_all # residuals of the fit
    fit = [missing_val]*n_all       # fitted regression line
    valid_count = 0          # total number of non-missing values used
    r_sum_sqr_x = 0.0        # sum square error between each value and the median
                             # of all the data
    r_sum_sqr_e = 0.0        # sum square error of residuals for each segment
                             # regression
                             # r_slope - slope of linear regression line
                             # r_t - slope error
    for i in range(n_all):
        if all_data[i] != missing_val:
            valid_count = valid_count + 1
            if valid_count < n_left:
                y_int = left_y_int
            else:
                y_int = right_y_int
            residuals[i] = (y_int + slope_est*(i+1)) - all_data[i]
            fit[i] = y_int + slope_est*(i+1)
            r_sum_sqr_e = r_sum_sqr_e + residuals[i]**2
            r_sum_sqr_x = r_sum_sqr_x + (float(i+1) - data_med)**2
    
    r_se_sqr = r_sum_sqr_e / (valid_count - 2)
    r_sb = sqrt(r_se_sqr / r_sum_sqr_x)
    r_t = slope_est / r_sb # unused?
    
    # Finally, we have fit same-sloped models to the left and right
    # segments. We can carry on the model test.
    r_mu = (left_y_int, right_y_int)
    r_alpha = slope_est
    sse_resid = r_sum_sqr_e
    
    # 3) Now it looks like we compute residuals for all our data using
    # the models attained above. I'm not sure why the PHA does this - it seems
    # to duplicate the error computed in the final step above in (2). Will
    # retain for now (TODO: rewrite at some point)
    r_residuals = [missing_val]*num_obs
    rssx = 0.0 # sum square error on the residuals
    rsse = [0.0, 0.0] # sum square error between each value and the mean for all
                      # the data.
    for k in range(num_obs):
        if all_data[k] != missing_val:
            if k < (bp_index+1):
                ind = 0
            else:
                ind = 1
            r_residuals[k] = all_data[k] - r_mu[ind] - r_alpha*(k+1)
            rsse[ind] = rsse[ind] + r_residuals[k]**2
            rssx = rssx + (float(k+1) - all_mean)**2
    
    # 4) We now have the squared error and are basically done!
    r_sum_sqr_tot = sum(rsse) # sum error on residuals for left and right
                              # segments
    
    # At this point, we have a few things - 
    #    r_mu - the y_intercepts of each segment for the model we fit
    #    r_alpha - the slope of the regression, based on the slopes generated
    #            in each left and right segment
    #    r_sum_sqr_tot - sum square total error of the residuals for each
    #                    segment.
    #
    # We now print out the info and compute critical values, BIC
    count = n_all
    sse_total = r_sum_sqr_tot # just computed
    
    # Compute F-statistic as the test statistic. The computation of this
    # statistic necessitates us to recall the error on the residuals we got
    # in KTHSLR1, when we assumed that there *was not* a change in slope at 
    # the suspected breakpoint
    f_val = ((sse_resid-sse_total)/1.)/(sse_total/(count-3))
    f_crit = lookup_critical(count-3, "f1")
    b_tot, b_err, b_pen = bayes(count, sse_total, 4)
    # amplitude change estimate
    y1 = r_mu[0] + r_alpha * range_all[bp_index+1] # end of left seg
    y2 = r_mu[1] + r_alpha * range_all[bp_index+2] # beg of right seg
    amp_est = y1-y2
    # k, we have finished this god-awful changepoint
    # output string
    head = "%7s %6.2f %7.2f %7.2f" % (cmodel, b_tot, b_err, b_pen)
    stats =  " %7.2f %7.2f %7.3f ------- %7.2f %7.2f" % (r_mu[0], r_mu[1],
                                                         r_alpha, f_val, f_crit)
    tail = "% 7.2f %5d %4d" % (amp_est, n_left, n_right)
    print (head+stats+tail)
    
    ## add to results dictionary
    changepoint_dict[cmodel] = dict(bic=b_tot, 
                                    mu=r_mu,
                                    alpha=[r_alpha, r_alpha],
                                    sse_bic=sqrt(sse_total/n_vals), # total number of non-missing
                                                                    # values, from KTHSLR0
                                    seg_lens=[n_left, n_right],
                                    test_stat=f_val,
                                    crit_val=f_crit,
                                    amp_est=amp_est)
    
    ############################################################################   
    ## TEST 4 - KTHTPR2, Amplitude shift, non-equally sloped segments
    ##    Method - Use Kendall-Theill method to compute a linear regression for
    ##            each of the left and right segments of the data, but do not
    ##            assume that they must have equal slopes. 
    cmodel = "KTHTPR2"
    
    all_data = y
    num_obs = len(y)
    
    ## NOTE - Using kth_line() to regress the left/right segments is already done
    ## in the test for KTHTPR0. Should be able to re-use the results from there
    ## in this test.
    left_x = x[:bp_index+1]
    left_y = y[:bp_index+1]
    right_x = x[bp_index+1:]
    right_y = y[bp_index+1:]
    
    kthl_left = kth_line(left_x, left_y, missing_val)
    kthl_right = kth_line(right_x, right_y, missing_val)
    ## END NOTE
    
    ## Estimate step change amplitude
    y1 = kthl_left.y_int + kthl_left.slope*left_x[-1]
    y2 = kthl_right.y_int + kthl_right.slope*right_x[0]
    amp_est = y1 - y2
    
    count = kthl_left.nval + kthl_right.nval
    #count = len(get_valid_data(left_data, missing_val))+len(get_valid_data(right_data, missing_val))
    
    sse_total = kthl_left.sseslope + kthl_right.sseslope
    print kthl_left.sseslope, kthl_right.sseslope, sse_resid
    
    # Compute test statistic, f2
    # Recall the residuals computed from KTHSLR1, where we assumed constant
    # slope throughout the data segment (sse_resid)
    f2_val = ((sse_resid-sse_total)/2.)/(sse_total/(count-4))
    f2_crit = lookup_critical(count-4, "f2")
    b_tot, b_err, b_pen = bayes(count, sse_total, 5)
    # output string
    head = "%7s %6.2f %7.2f %7.2f" % (cmodel, b_tot, b_err, b_pen)
    stats =  " %7.2f %7.2f %7.3f %7.3f %7.2f %7.2f" % (kthl_left.y_int, kthl_right.y_int,
                                                       kthl_left.slope, kthl_right.slope,
                                                       f2_val, f2_crit)
    tail = "% 7.2f %5d %4d" % (amp_est, n_left, n_right)
    print (head+stats+tail)       
    
    ## add to results dictionary
    changepoint_dict[cmodel] = dict(bic=b_tot, 
                                    mu=[kthl_left.y_int, kthl_right.y_int],
                                    alpha=[kthl_left.slope, kthl_right.slope],
                                    sse_bic=sqrt(sse_total/n_vals), # total number of non-missing
                                                                    # values, from KTHSLR0
                                    seg_lens=[n_left, n_right],
                                    test_stat=f2_val,
                                    crit_val=f2_crit,
                                    amp_est=amp_est)
       
    ############################################################################   
    ## TEST 5 - KTHTPR3, Amplitude shift, flat-to-sloped segments
    ##    Method - Assume - using the regressions fit to each left/right segment
    ##            previously - that the segment on the left has slope = 0.0.
    cmodel = "KTHTPR3"
    
    # We will re-use the kthl_left/right regressions fit in the last model
    y1 = kthl_left.y_med
    y2 = kthl_right.y_int + kthl_right.slope*right_x[0]
    amp_est = y1 - y2
    
    count = kthl_left.nval + kthl_right.nval
    #count = len(get_valid_data(left_data, MISS))+len(get_valid_data(right_data, MISS))
    
    sse_total = kthl_left.sseflat + kthl_right.sseslope
    # Compute test statistic, f1
    # Recall the residuals computed from KTHSLR1, where we assumed constant
    # slope throughout the data segment (sse_resid)
    f_val =  ((sse_resid-sse_total)/1.0)/(sse_total/(count-3))
    f_crit = lookup_critical(count, "f1")
    b_tot, b_err, b_pen = bayes(count, sse_total, 4)
    # output string
    head = "%7s %6.2f %7.2f %7.2f" % (cmodel, b_tot, b_err, b_pen)
    stats =  " %7.2f %7.2f ------- %7.3f %7.2f %7.2f" % (kthl_left.y_med, kthl_right.y_int,
                                                        kthl_right.slope, f_val, f_crit)
    tail = "% 7.2f %5d %4d" % (amp_est, n_left, n_right)
    print (head+stats+tail)
    
    ## add to results dictionary
    changepoint_dict[cmodel] = dict(bic=b_tot, 
                                    mu=[kthl_left.y_med, kthl_right.y_int],
                                    alpha=[0.0, kthl_right.slope],
                                    sse_bic=sqrt(sse_total/n_vals), # total number of non-missing
                                                                    # values, from KTHSLR0
                                    seg_lens=[n_left, n_right],
                                    test_stat=f_val,
                                    crit_val=f_crit,
                                    amp_est=amp_est)
    
    ############################################################################   
    ## TEST 6 - KTHTPR4, Amplitude shift, sloped-to-flat segments
    ##    Method - Assume - using the regressions fit to each left/right segment
    ##            previously - that the segment on the right has slope = 0.0.
    cmodel = "KTHTPR4"
    
    # We will re-use the kthl_left/right regressions fit in KTHTPR2
    y1 = kthl_left.y_int + kthl_left.slope*left_x[-1]
    y2 = kthl_right.y_med
    amp_est = y1 - y2
    
    count = kthl_left.nval + kthl_right.nval    
    #count = len(get_valid_data(left_data, MISS))+len(get_valid_data(right_data, MISS))
    
    # Compute test statistic, f1
    # Recall the residuals computed from KTHSLR1, where we assumed constant
    # slope throughout the data segment (sse_resid)
    sse_total = kthl_left.sseslope + kthl_right.sseflat
    f_val = ((sse_resid-sse_total)/1.0)/(sse_total/(count-3))
    f_crit = lookup_critical(count, "f1")
    b_tot, b_err, b_pen = bayes(count, sse_total, 4)
    # output string
    head = "%7s %6.2f %7.2f %7.2f" % (cmodel, b_tot, b_err, b_pen)
    stats =  " %7.2f %7.2f %7.3f ------- %7.2f %7.2f" % (kthl_left.y_int, kthl_right.y_med,
                                                         kthl_left.slope, f_val, f_crit)
    tail = "% 7.2f %5d %4d" % (amp_est, n_left, n_right)
    print (head+stats+tail)
    
    ## add to results dictionary
    changepoint_dict[cmodel] = dict(bic=b_tot, 
                                    mu=[kthl_left.y_int, kthl_right.y_med],
                                    alpha=[kthl_left.slope, 0.0],
                                    sse_bic=sqrt(sse_total/n_vals), # total number of non-missing
                                                                    # values, from KTHSLR0
                                    seg_lens=[n_left, n_right],
                                    test_stat=f_val,
                                    crit_val=f_crit,
                                    amp_est=amp_est)
    
    ################################################################################ 
    ## We've finished the basic TPR tests. The PHA has loop-back criteria
    ## which could cause the code to double back and compute these test
    ## statistics over again or with subtle variations. I'll come back to
    ## that...
    ##
    ## For now, let's look at the test results and pick the most likely
    ## changepoint class based on the minimum BIC value found
    results = [(cmodel, changepoint_dict[cmodel]['bic']) for cmodel in changepoint_dict]
    results = sorted(results, key=itemgetter(1))
    cmodel, q = results[0]
        
    output = changepoint_dict[cmodel]
    iqtype = minbic_models.index(cmodel) + 1
    sseq=output['sse_bic']
    rmu=output['mu']
    slpq=output['alpha']
    qstat=output['test_stat']
    qcrit=output['crit_val']
    qoff=output['amp_est']
    mkfinq=output['seg_lens']
    
    print ("Post: - %2d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %5d %5d" %
           (iqtype, q, sseq, rmu[0], rmu[1], slpq[0], slpq[1], qstat, qcrit, qoff,
            mkfinq[0], mkfinq[1]) )
    
    ## Store information about this breakpoint and its crucial stats
    return dict(cmodel=cmodel, iqtype=iqtype, asigx=qoff, azscr=qoff/sseq, 
                rslp=[kthl_left.slope, kthl_right.slope])
    
        
    