#!/usr/bin/env python
#
# Copyright (C) 2011 Daniel Rothenberg.
# See Google Code project page for license, 
# http://code.google.com/p/ccf-homogenization/

"""Port of functions necessary for identifying and classifying changepoint
model types during splitmerge process

GLOSSARY OF ABBREVIATIONS -

breakpoint/changepoint - Same meaning; an index in a dataset where it is hypothesized
    that either the mean in the data or the linear trend in the data changes.
BIC - Bayes Information Criterion; test statistic for evaluating test statistics
KTH - Kendall-Theill robust line fit method 
    (http://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator)
TPR - two-phase regression; the process of fitting separate linear regressions
    to a set of data on each side of a suspected discontinuity
SLR - Straight-line regression; the normal process of fitting a linear
    regression to a set of data

"""
__docformat__ = "restructuredtext"

# http://docs.python.org/library/math.html
from math import log10, sqrt
# http://docs.python.org/library/operator.html
import operator

import numpy as np

# ccf-homogenization imports
from util import compute_mean, get_valid_data, compute_std, median

class Stats(dict):
    """Base object for returning results of statistical tests.
    
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
    def valid_pair(pair, missing_val=-9999):
        return ( (pair[0] != missing_val))
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
         
    # calculate slope and y intercept   
    numer = 0.0
    denom = 0.0
    for (x_val, y_val) in good_data:
        numer = numer + (x_val-x_mean)*(y_val-y_mean)
        denom = denom + (x_val-x_mean)*(x_val-x_mean)
    
    slope = numer/denom
    y_int = y_mean - slope*x_mean
    
    # calculate residuals
    sseflat = 0.0
    sseslope = 0.0
    for (x_val, y_val) in good_data:
        sseflat = sseflat + (y_mean-y_val)**2
        resid = (y_int + slope*x_val) - y_val
        sseslope = sseslope + resid**2
    
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
    def valid_pair(pair, missing_val=-9999):
        """Given a 2-tuple, returns False if either of the values are equal
        to the supplied missing_val.
        """
        return ( (pair[0] != missing_val) and 
                 (pair[1] != missing_val))
    
    pairs = zip(x, y)
    # Filter out data pairs where either x_i or y_i equal missing_val
    good_data = [pair for pair in pairs if valid_pair(pair)]
    #good_x, good_y = zip(*good_data)
    good_x, good_y = map(list, zip(*good_data))
    nval = len(good_data)
        
    # calculate paired slopes, find median value
    if nval > 0:
        nslp = 0
        slopes = []
        for i in range(nval-1):
            for j in range(i, nval):
                if good_x[i] != good_x[j]:
                    nslp += 1
                    slopes.append( (good_y[j]-good_y[i])/(good_x[j]-good_x[i]) )   
        slope = np.median(slopes)
        
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
#    
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

def minbic(x, y, breakpoint, missing_val=-9999, models=None):
    """Performs a changepoint-detection analysis to classify a suspected
    changepoint as one of a given list of model types.
    
    This method performs several statistical tests in an attempt to fit a
    model to the suspected changepoint in the data passed to it.
    
    Upon completion of the changepoint tests, this method will identify which
    test minimizes the Bayesian Information Criterion, based on the error
    computed for each test. It will then return that changepoint's statistics
    for later use.
    
    :Params x, y:
        The predictor and predictand dataset where there is a suspected
        changepoint. Must be equal length lists of floats, but may contain
        missing_val in place of missing values.
    :Param breakpoint:
        The index of the suspected changepoint in the x and y datasets. Must be
        in range(len(x))
    :Param missing_val:
        (optional) A float or int representing the placeholder for missing 
        values in the datasets; default is -9999.
    :Param models:
        The list of potential model types and analysis functions. Should be 
        passed as a list of 2-tuples. The first element in each tuple should be
        the model name as a string; the second element should be the function
        to call to analyze this model. For example, a valid models passed to
        this function might be:
            [('my_simple_test', test_function), ...]
    :Return:
        A dictionary containing the results of this analysis with the following
        fields:
        :Param model_name:
            The name of the changepoint model, as a string
        :Param iqtype:
            The type of changepoint model, as an integer. Corresponds to the
            order in the 'models' list that was passed to this function.
        :Param offset:
            The computed amplitude change at the breakpoint.
        :Param offset_z:
            The z-score of the computed amplitude change at the breakpoint;
            the amplitude normalized by the error from the model fit 
            at the breakpoint.
        :Param slopes:
            A two-element tuple containing the slope of the regression fit to
            the left and right segments of the data.
            
    ...
    
    :Raises AssertionError:
        If the length of x and y are not equal.
            
    """    
    ## Figure out what changepoint models to test for. A default is defined
    ## here if no models were specified by the user.
    if not models:
        models = [
                # Homogeneous linear types -
                   #('KTHSLR0', kthslr0) # simple flat line
                    ('KTHSLR1', kthslr1), # sloped straight line
                 # Two phase regression types - 
                    ('KTHTPR0', kthtpr0), # amplitude only, 0 sloped segments (step change)
                    ('KTHTPR1', kthtpr1), # amplitude shift with equal sloped segments
                    ('KTHTPR2', kthtpr2), # amplitude shift with non-equal sloped segments
                    ('KTHTPR3', kthtpr3), # amplitude shift with flat-to-slope segments
                    ('KTHTPR4', kthtpr4), # amplitude shift with slope-to-flat segments  
                 ]
    model_order, _ = zip(*models)
        
    # Sanity check that we have a good batch of data
    assert len(y) == len(x)

    ## All of the changepoint model tests which will be performed here print
    ## a one-line debug summary in a formatted table; this is the header row for
    ## that table.
    left_header = " QTYP     QVAL    QRSE     QPF     MU1     MU2  ALPHA1"
    right_header = "  ALPHA2   MSTAT   MCRIT    MOFF KNT1 KNT2"
    #print (left_header+right_header)
    
    # The following container class will help us record computed values (like
    # errors for different models, etc) for re-use in the following tests.
    class Values:
        """Container class for recording and passing along values computed by
        each of the changepoint detection functions.
                
        """
        pass
    vals = Values()
    
    # For storing the results of these tests.
    changepoint_dict = {}
    
    ## Loop over the models we're using to classify this breakpoint and perform
    ## the specified analysis
    for (model_name, model_test) in models:
        result = model_test(x, y, breakpoint, missing_val, vals)
        changepoint_dict[model_name] = result
       
    ############################################################################ 
    ## Identify the best-fitting model computed above by finding which one has
    ## the lowest BIC value
    results = [(model, changepoint_dict[model]['bic']) for (model,_) in models]
    results = sorted(results, key=operator.itemgetter(1))
        
    cmodel, bic = results[0]
    output = changepoint_dict[cmodel]
    
    ############################################################################ 
    ## Debug printing for comparison to Fortran PHA        
    iqtype = model_order.index(cmodel) + 2
    sse_bic=output['sse_bic']
    mu=output['mu']
    slp=output['alpha']
    test_stat=output['test_stat']
    crit_val=output['crit_val']
    offset=output['amp_est']
    seg_lens=output['seg_lens']
    
    #print ("Post: - %2d %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %5d %5d" %
    #       (iqtype, bic, sse_bic, mu[0], mu[1], slp[0], slp[1], test_stat, crit_val, offset,
    #        seg_lens[0], seg_lens[1]) )
    ############################################################################ 
    
    ## Return the details of the best-fitting model
    # The slopes from each left/right segment
    kthl_left = vals.left_kth_regress
    kthl_right = vals.right_kth_regress
    slopes = [kthl_left.slope, kthl_right.slope] 
    # The amplitude change. We need to use the simple one from 'KTHTPR0' if a
    # more complex two-phase regression was the best fit - we want the step 
    # change in the median values of the data on each side. This is a more 
    # robust estimate of the center of the data than the mean is and let's us
    # see how big the jump is before and after the breakpoint.
    if 'TPR' in cmodel:
        offset = kthl_left.y_med - kthl_right.y_med
    else:
        offset = output['amp_est']
    # Z-score for offset, normalized about the sum square error of the residuals
    # of the model fit.
    offset_z = offset/sse_bic
    
    ## Store information about this breakpoint and its crucial stats
    return dict(cmodel=cmodel, iqtype=iqtype, offset=offset, offset_z=offset_z, 
                slopes=slopes)

def kthslr1(x, y, breakpoint, missing_val, vals):
    """KTHSLR1, Sloped straight line test
    
    Utilizes the Kendall-Theill method to compute a linear regression on all
    the given data. With this regression, it computes the error for this test
    as the sum square error of modeled residuals.
    
    :Params x, y:
        The predictor and predictand datasets as lists of floats. The length of
        x and y must be the same size, but the values may contain an arbitrary
        number of missing values as specified by missing_val.
    :Param breakpoint:
        The index of the suspected changepoint in the data. Must be
        in range(len(x))
    :Param missing_val:
        A float or int representing the placeholder value for missing
        values in the dataset
    :Param vals:
        A 'Values' object for storing important intermediate values from the 
        computations during this test. This test does not require vals to 
        contain any previously computed values.
    :Returns:
        result, a dictionary with the following fields:
            :Param cmodel:
                The name of the test ('KTHSLR1').
            :Param bic: 
                The Bayes Information Criterion statistic for this analysis.
            :Param mu: 
                A 2-value tuple with the y-intercepts of the regression here. 
                Since the y-intercept computed here is for the regression fit 
                to the entire dataset, mu will be [computed_intercept, computed_intercept].
            :Param alpha: 
                A 2-value tuple with the computed slopes of the regressions 
                here. Since the regression fit here is for the entire data set,
                alpha will be [computed_slope, computed_slope]
            :Param sse_bic: 
                The sum square error on the computed BIC values.
            :Param seg_lens: 
                A 2-value tuple with the length of the segments in the 
                regressions here. Since we consider the whole dataset here, will
                be [total_length, 0].
            :Param test_stat: 
                The test statistic value computed here. 0.0 by default since no
                test is being performed.
            :Param crit_val: 
                The critical value for the test statistic at alpha=0.05. 0.0 by
                default since no test is being performed.
            :Param amp_est: 
                An estimate of the changepoint amplitude. 0.0 by default since
                we are assuming there is no step change in the data.
        vals, the original 'Values' object provided to this method with the
            following modifications:
            added vals.sse_resid - the sum square error on the residuals from this
                regression
            added vals.n_vals - the number of valid (x, y) pairs found found in
                the regression here
            added vals.all_kth_regress - the regression results when the kth
                method was used to fit a linear regression to all the (x, y) data.
                
    """
    cmodel = "KTHSLR1"
    kth_regression = kth_line(x, y, missing_val)
    
    num_obs =  len(x) # total number of observations
    n_vals = kth_regression.nval # number of non-missing observations
    slope = kth_regression.slope
    y_int = kth_regression.y_int
    sse_resid = kth_regression.sseslope
    amp_est = 0.0 # Test assumes no step change
    
    # Bayes Info Criterion - total, due to error in variance, and penalty
    bic, bic_error, bic_dof_penalty = bayes(num_obs, sse_resid, 2)
    # output string
    head = "%7s %6.2f %7.2f %7.2f" % (cmodel, bic, bic_error, bic_dof_penalty)
    stats =  " %7.2f ------- %7.3f ------- ------- -------" % (y_int, slope)
    tail = " %7.2f %5d ----" % (amp_est, n_vals)
    #print (head+stats+tail)
    
    ## Record the results of this test
    result = dict(cmodel=cmodel,
                  bic=bic, 
                  mu=[y_int, y_int],
                  alpha=[slope, slope], 
                  sse_bic=(sse_resid/n_vals),
                  seg_lens=[n_vals, 0],
                  test_stat=0.0,
                  crit_val=0.0,
                  amp_est=amp_est)
    
    ## Add the important intermediate values from this analysis to our
    ## 'vals' object
    vals.sse_resid = sse_resid
    vals.n_vals = n_vals
    vals.all_kth_regress = kth_regression
    
    return result, vals

def kthtpr0(x, y, breakpoint, missing_val, vals):
    """KTHTPR0, Amplitude change, no slope change
    
    Utilizes the Kendall-Theill method to compute a linear regression on the 
    data segments on each side of the suspected breakpoint. Then, applies a
    Student's t-test to see whether the mean value for the two segments differs.
    
    :Params x, y:
        The predictor and predictand datasets as lists of floats. The length of
        x and y must be the same size, but the values may contain an arbitrary
        number of missing values as specified by missing_val.
    :Param breakpoint:
        The index of the suspected changepoint in the data. Must be
        in range(len(x))
    :Param missing_val:
        A float or int representing the placeholder value for missing
        values in the dataset
    :Param vals:
        A 'Values' object for storing important intermediate values from the 
        computations during this test. This test does not require vals to 
        contain any previously computed values.
    :Returns:
        result, a dictionary with the following fields:
            :Param cmodel:
                The name of the test ('KTHTPR0')
            :Param bic: 
                The Bayes Information Criterion statistic for this analysis.
            :Param mu: 
                A 2-value tuple with the y-intercepts of the regression here.
                Should be [left segment y-intercept, right segment y-intercept].
            :Param alpha: 
                A 2-value tuple with the computed slopes of the regressions 
                here. Should be [left segment slope, right segment slope].
            :Param sse_bic: 
                The sum square error on the computed BIC value.
            :Param seg_lens: 
                A 2-value tuple with the length of the segments in the 
                regressions here. Should be [length of left segment, length of
                right segment].
            :Param test_stat: 
                The test statistic value computed here.
            :Param crit_val: 
                The critical value for the test statistic at alpha=0.05.
            :Param amp_est: 
                An estimate of the changepoint amplitude based on the medians
                of the left and right segments.
        vals, the original 'Values' object provided to this method with the
            following modifications:
            added vals.left_kth_regress - the regression fit on the left 
                segment of data
            added vals.n_left - the number of valid (x, y) data pairs in the 
                left segment of data
            added vals.right_kth_regress - the regression fit on the right 
                segment of data
            added vals.n_right - the number of valid (x, y) data pairs in the 
                right segment of data.
            added vals.tpr_offset - the step change difference in the mean from
                the left to the right segments of data
                
    """
    cmodel = "KTHTPR0"
    
    # Divide the segment before/after the suspect breakpoint
    left_x = x[:breakpoint+1]
    left_y = y[:breakpoint+1]
    right_x = x[breakpoint+1:]
    right_y = y[breakpoint+1:]
    
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
    bic, bic_error, bic_dof_penalty = bayes(n_total, sse_sum, 3)
    # output string
    head = "%7s %6.2f %7.2f %7.2f" % (cmodel, bic, bic_error, bic_dof_penalty)
    stats =  " %7.2f %7.2f ------- ------- %7.2f %7.2f" % (left_y_med, right_y_med,
                                                          t_val, t_crit)
    tail = "% 7.2f %5d %4d" % (amp_est, n_left, n_right)
    #print (head+stats+tail)
    
    ## Record the results of this test
    result = dict(cmodel=cmodel,
                  bic=bic, 
                  mu=[left_y_med, right_y_med],
                  alpha=[0.0, 0.0],
                  sse_bic=sqrt(sse_sum/n_total),
                  seg_lens=[n_left, n_right],
                  test_stat=t_val,
                  crit_val=t_crit,
                  amp_est=amp_est)
    
    ## Add the important intermediate values from this analysis to our
    ## 'vals' object
    vals.left_kth_regress = kthl_left
    vals.n_left = n_left
    vals.right_kth_regress = kthl_right
    vals.n_right = n_right
    vals.tpr_offset = amp_est
    
    return result

def kthtpr1(x, y, breakpoint, missing_val, vals):
    """KTHTPR1, Amplitude change, equal slope segments
    
    Utilizes the Kendall-Theill method to compute a linear regression on the 
    data segments on each side of the suspected breakpoint. Then, assuming the
    slope on each segment is equal, computes errors from the modeled residuals
    on each segment.
    
    :Params x, y:
        The predictor and predictand datasets as lists of floats. The length of
        x and y must be the same size, but the values may contain an arbitrary
        number of missing values as specified by missing_val.
    :Param breakpoint: 
        The index of the suspected changepoint in the data. Must be
        in range(len(x))
    :Param missing_val:
        A float or int representing the placeholder value for missing
        values in the dataset
    :Param vals:
        A 'Values' object for storing important intermediate values from the 
        computations during this test. This test requires vals to already 
        contain to values:
            vals.sse_resid - the sum square error on the residuals computed by
                fitting a linear regression on the entire dataset
            vals.n_vals - the number of valid (x, y) pairs given by the user
    :Returns:
        result, a dictionary with the following fields:
            :Param cmodel:
                The name of the test ('KTHTPR1').
            :Param bic: 
                The Bayes Information Criterion statistic for this analysis.
            :Param mu: 
                A 2-value tuple with the y-intercepts of the regression fit to
                the segment on each side of the breakpoint.
            :Param alpha: 
                A 2-value tuple with the computed slopes of the regressions 
                here fit on each side of the breakpoint.
            :Param sse_bic: 
                The sum square error on the computed BIC value.
            :Param seg_lens: 
                A 2-value tuple with the length of the segments in the 
                regressions here on each side of the breakpoint.
            :Param test_stat: 
                The test statistic value computed here.
            :Param crit_val: 
                The critical value for the test statistic at alpha=0.05.
            :Param amp_est: 
                An estimate of the changepoint amplitude.
        vals, the original 'Values' object provided to this method with no
            modifications.
                
    """
    cmodel = "KTHTPR1"
    
    num_obs = len(y)
    left_data = y[:breakpoint+1]
    right_data = y[breakpoint+1:]
    
    ## 1) Compute the mean for *all* of the data
    all_valid_data = get_valid_data(y, missing_val)
    all_mean = compute_mean(all_valid_data, valid=True)
    
    ## 2) use Kendall-Theill method with single slope
    # This method is slightly different than the kth_line() method above.
    # First, we get only the valid data, and we pair it with the natural ordering
    # of the data, i.e. 1, 2, 3...
    ## BUG: The computation for range_all here will result in a the predictor
    ##    dataset ultimately being equal to the order of the valid data. To
    ##    illustrate this problem, consider the x data [1,2,3,4,5,6] and the y
    ##    data [10, -9999, -9999, -9999, 4]. With this method, we will ultimately
    ##    compute paired slopes with the points (1, 10) and (2, 4) when we really
    ##    need to use the points (1, 10) and (6, 4). The difference is signficant,
    ##    but it is an easy change - merely reproduce the valid_x and valid_y
    ##    code from kth_line() above.
    n_all = len(all_valid_data)
    range_all = range(1, n_all+1)
    
    valid_left = get_valid_data(left_data, missing_val)
    valid_right = get_valid_data(right_data, missing_val)
    
    n_left, n_right = len(valid_left), len(valid_right)
    range_left = range(1, n_left+1)
    range_right = range(n_left+1, n_all+1)
    
    # Second, generate paired slopes for the first segment.
    slopes = []
    for i in range(n_left-1):
        for j in range(i, n_left):
            if range_left[j] != range_left[i]:
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
                slopes.append( (valid_right[j]-valid_right[i])/
                               (range_right[j]-range_right[i]) )
    nslp = len(slopes)
    #slopes = sorted(slopes)

    # BUG: This is the median-index finding method used in the original
    #    PHA. It's different from what is used in kth_line() (which is
    #    the correct way to find the median in a list of data - average the
    #    middle two values if the list length is even).
    slope_est = median(slopes)
        
    #print "slope, ic, imed: %7.2f %5d %5d" % (slope_est, nslp, imed)
    
    # Fifth, compute the first segment intercept, y-median - slope*x-median
    # BUG: Another instance where we do not handle the case of an even number of
    #    values when computing the median
    imed = (n_left - 1)/2
    if (n_left%2)==1: imed = imed+1
    range_med = range_left[imed]
    valid_left = sorted(valid_left)
    data_med = valid_left[imed]
    left_y_int = data_med-slope_est*range_med
    #print "Seg1 - Xmed, Ymed, slope, Yint: %7.2f %7.2f %7.2f %7.3f" % (range_med, 
    #                                                                   data_med,
    #                                                                   slope_est,
    #                                                                   left_y_int)
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
    #print "Seg2 - Xmed, Ymed, slope, Yint: %7.2f %7.2f %7.2f %7.3f" % (range_med, 
    #                                                                   data_med,
    #                                                                   slope_est,
    #                                                                   right_y_int)
    
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
        if all_valid_data[i] != missing_val:
            valid_count = valid_count + 1
            if valid_count < n_left:
                y_int = left_y_int
            else:
                y_int = right_y_int
            residuals[i] = (y_int + slope_est*(i+1)) - all_valid_data[i]
            fit[i] = y_int + slope_est*(i+1)
            r_sum_sqr_e += r_sum_sqr_e + residuals[i]**2
            r_sum_sqr_x += r_sum_sqr_x + (float(i+1) - data_med)**2
    
    r_se_sqr = r_sum_sqr_e / (valid_count - 2)
    r_sb = sqrt(r_se_sqr / r_sum_sqr_x)
    r_t = slope_est / r_sb # unused?
    
    # Finally, we have fit same-sloped models to the left and right
    # segments. We can carry on the model test.
    r_mu = (left_y_int, right_y_int)
    r_alpha = slope_est
    sse_seg_resid = r_sum_sqr_e
    
    # 3) Now it looks like we compute residuals for all our data using
    # the models attained above. I'm not sure why the PHA does this - it seems
    # to duplicate the error computed in the final step above in (2). Will
    # retain for now (TODO: rewrite at some point)
    rssx = 0.0 # sum square error on the residuals
    rsse = [0.0, 0.0] # sum square error between each value and the mean for all
                      # the data.
    for k in range(num_obs):
        if all_valid_data[k] != missing_val:
            if k <= breakpoint:
                ind = 0
            else:
                ind = 1
            residual = all_valid_data[k] - r_mu[ind] - r_alpha*(k+1)
            rsse[ind] = rsse[ind] + residual**2
            rssx += (float(k+1) - all_mean)**2
    
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
    f_val = ((vals.sse_resid-sse_total)/1.)/(sse_total/(count-3))
    f_crit = lookup_critical(count-3, "f1")
    bic, bic_error, bic_dof_penalty = bayes(count, sse_total, 4)
    # amplitude change estimate
    #print r_mu, r_alpha, breakpoint
    #print len(range_all)
    #y1 = r_mu[0] + r_alpha * range_all[breakpoint+1] # end of left seg
    #y2 = r_mu[1] + r_alpha * range_all[breakpoint+2] # beg of right seg
    y1 = r_mu[0] + r_alpha * range_left[-1]
    y2 = r_mu[1] + r_alpha * range_right[0]
    amp_est = y1-y2
    # k, we have finished this god-awful changepoint
    # output string
    head = "%7s %6.2f %7.2f %7.2f" % (cmodel, bic, bic_error, bic_dof_penalty)
    stats =  " %7.2f %7.2f %7.3f ------- %7.2f %7.2f" % (r_mu[0], r_mu[1],
                                                         r_alpha, f_val, f_crit)
    tail = "% 7.2f %5d %4d" % (amp_est, n_left, n_right)
    #print (head+stats+tail)
    
    ## Record results of this test
    result = dict(cmodel=cmodel,
                  bic=bic, 
                  mu=r_mu,
                  alpha=[r_alpha, r_alpha],
                  sse_bic=sqrt(sse_total/vals.n_vals), # total number of non-missing values, from KTHSLR0
                  seg_lens=[n_left, n_right],
                  test_stat=f_val,
                  crit_val=f_crit,
                  amp_est=amp_est)
    
    return result

def kthtpr2(x, y, breakpoint, missing_val, vals):
    """KTHTPR2, Amplitude shift, non-equally sloped segments
    
    Utilizes the Kendall-Theill method to compute a linear regression on the 
    data segments on each side of the suspected breakpoint. Then, uses these 
    two regressions to compute errors on the residuals from each data segment.
    
    :Params x, y:
        The predictor and predictand datasets as lists of floats. The length of
        x and y must be the same size, but the values may contain an arbitrary
        number of missing values as specified by missing_val.
    :Param breakpoint:
        The index of the suspected changepoint in the data. Must be
        in range(len(x))
    :Param missing_val:
        A float or int representing the placeholder value for missing
        values in the dataset
    :Param vals:
        A 'Values' object for storing important intermediate values from the 
        computations during this test. This test requires vals to already 
        contain to values:
            vals.sse_resid - the sum square error on the residuals computed by
                fitting a linear regression on the entire dataset
            vals.n_vals - the number of valid (x, y) pairs given by the user
            vals.left_kth_regress - the regression on the left data segment
            vals.n_left - the number of valid (x, y) pairs in the left segment
            vals.right_kth_regress - the regression on the right data segment
            vals.n_right - the number of valid (x, y) pairs in the right segment
    :Returns:
        result, a dictionary with the following fields:
            :Param cmodel:
                The name of the test ('KTHTPR2').
            :Param bic: 
                The Bayes Information Criterion statistic for this analysis.
            :Param mu: 
                A 2-value tuple with the y-intercepts of the regression here.
            :Param alpha: 
                A 2-value tuple with the computed slopes of the regressions .
                here.
            :Param sse_bic: 
                The sum square error on the computed BIC value.
            :Param seg_lens: 
                A 2-value tuple with the length of the segments in the 
                regressions here.
            :Param test_stat: 
                The test statistic value computed here.
            :Param crit_val: 
                The critical value for the test statistic at alpha=0.05.
            :Param amp_est: 
                An estimate of the changepoint amplitude.
        vals, the original 'Values' object provided to this method with no
            modifications.
                
    """
    cmodel = "KTHTPR2"
        
    ## NOTE - Using kth_line() to regress the left/right segments is already done
    ## in the test for KTHTPR0. Should be able to re-use the results from there
    ## in this test.
    left_x = x[:breakpoint+1]
    right_x = x[breakpoint+1:]
    
    #kthl_left = kth_line(left_x, left_y, missing_val)
    #kthl_right = kth_line(right_x, right_y, missing_val)
    kthl_left = vals.left_kth_regress
    kthl_right = vals.right_kth_regress
    ## END NOTE
    
    ## Estimate step change amplitude
    y1 = kthl_left.y_int + kthl_left.slope*left_x[-1]
    y2 = kthl_right.y_int + kthl_right.slope*right_x[0]
    amp_est = y1 - y2
    
    count = kthl_left.nval + kthl_right.nval
    
    sse_total = kthl_left.sseslope + kthl_right.sseslope
    #print kthl_left.sseslope, kthl_right.sseslope, vals.sse_resid
    
    # Compute test statistic, f2
    # Recall the residuals computed from KTHSLR1, where we assumed constant
    # slope throughout the data segment (sse_resid)
    f2_val = ((vals.sse_resid-sse_total)/2.)/(sse_total/(count-4))
    f2_crit = lookup_critical(count-4, "f2")
    bic, bic_error, bic_dof_penalty = bayes(count, sse_total, 5)
    # output string
    head = "%7s %6.2f %7.2f %7.2f" % (cmodel, bic, bic_error, bic_dof_penalty)
    stats =  " %7.2f %7.2f %7.3f %7.3f %7.2f %7.2f" % (kthl_left.y_int, kthl_right.y_int,
                                                       kthl_left.slope, kthl_right.slope,
                                                       f2_val, f2_crit)
    tail = "% 7.2f %5d %4d" % (amp_est, vals.n_left, vals.n_right)
    #print (head+stats+tail)       
    
    ## Record result of this test
    result = dict(cmodel=cmodel,
                  bic=bic, 
                  mu=[kthl_left.y_int, kthl_right.y_int],
                  alpha=[kthl_left.slope, kthl_right.slope],
                  sse_bic=sqrt(sse_total/vals.n_vals), # total number of non-missing values, from KTHSLR0
                  seg_lens=[vals.n_left, vals.n_right],
                  test_stat=f2_val,
                  crit_val=f2_crit,
                  amp_est=amp_est)
    
    return result

def kthtpr3(x, y, breakpoint, missing_val, vals):
    """KTHTPR3, Amplitude shift, flat-to-sloped segments
    
    Utilizes the Kendall-Theill method to compute a linear regression on the 
    data segments on each side of the suspected breakpoint. Then, uses these 
    two regressions to compute errors on the residuals from each data segment,
    while assuming that the slope in the left segment is actually 0.0.
    
    :Params x, y:
        The predictor and predictand datasets as lists of floats. The length of
        x and y must be the same size, but the values may contain an arbitrary
        number of missing values as specified by missing_val.
    :Param breakpoint:
        The index of the suspected changepoint in the data. Must be
        in range(len(x))
    :Param missing_val:
        A float or int representing the placeholder value for missing
        values in the dataset
    :Param vals:
        A 'Values' object for storing important intermediate values from the 
        computations during this test. This test requires vals to already 
        contain to values:
            vals.sse_resid - the sum square error on the residuals computed by
                fitting a linear regression on the entire dataset
            vals.n_vals - the number of valid (x, y) pairs given by the user
            vals.left_kth_regress - the regression on the left data segment
            vals.n_left - the number of valid (x, y) pairs in the left segment
            vals.right_kth_regress - the regression on the right data segment
            vals.n_right - the number of valid (x, y) pairs in the right segment
    :Returns:
        result, a dictionary with the following fields:
            :Param cmodel:
                The name of the test ('KTHTPR3').
            :Param bic: 
                The Bayes Information Criterion statistic for this analysis.
            :Param mu: 
                A 2-value tuple with the y-intercepts of the regression here.
            :Param alpha: 
                A 2-value tuple with the computed slopes of the regressions 
                here.
            :Param sse_bic: 
                The sum square error on the computed BIC value.
            :Param seg_lens: 
                A 2-value tuple with the length of the segments in the 
                regressions here.
            :Param test_stat: 
                The test statistic value computed here.
            :Param crit_val: 
                The critical value for the test statistic at alpha=0.05.
            :Param amp_est: 
                An estimate of the changepoint amplitude.
        vals, the original 'Values' object provided to this method with no
            modifications.
                
    """
    cmodel = "KTHTPR3"
    
    right_x = x[breakpoint+1:]
    
    ## Recall the regression fit on the left and right segments
    kthl_left = vals.left_kth_regress
    kthl_right = vals.right_kth_regress
    
    # We will re-use the kthl_left/right regressions fit in the last model
    y1 = kthl_left.y_med
    y2 = kthl_right.y_int + kthl_right.slope*right_x[0]
    amp_est = y1 - y2
    
    count = kthl_left.nval + kthl_right.nval
    
    sse_total = kthl_left.sseflat + kthl_right.sseslope
    # Compute test statistic, f1
    # Recall the residuals computed from KTHSLR1, where we assumed constant
    # slope throughout the data segment (sse_resid)
    f_val =  (vals.sse_resid-sse_total)/(sse_total/(count-3))
    f_crit = lookup_critical(count, "f1")
    bic, bic_error, bic_dof_penalty = bayes(count, sse_total, 4)
    # output string
    head = "%7s %6.2f %7.2f %7.2f" % (cmodel, bic, bic_error, bic_dof_penalty)
    stats =  " %7.2f %7.2f ------- %7.3f %7.2f %7.2f" % (kthl_left.y_med, kthl_right.y_int,
                                                        kthl_right.slope, f_val, f_crit)
    tail = "% 7.2f %5d %4d" % (amp_est, vals.n_left, vals.n_right)
    #print (head+stats+tail)
    
    ## Records results from this analysis
    result = dict(cmodel=cmodel,
                  bic=bic, 
                  mu=[kthl_left.y_med, kthl_right.y_int],
                  alpha=[0.0, kthl_right.slope],
                  sse_bic=sqrt(sse_total/vals.n_vals), # total number of non-missing values, from KTHSLR0
                  seg_lens=[vals.n_left, vals.n_right],
                  test_stat=f_val,
                  crit_val=f_crit,
                  amp_est=amp_est)
    
    return result

def kthtpr4(x, y, breakpoint, missing_val, vals):
    """KTHTPR4, Amplitude shift, sloped-to-flat segments
    
    Utilizes the Kendall-Theill method to compute a linear regression on the 
    data segments on each side of the suspected breakpoint. Then, uses these 
    two regressions to compute errors on the residuals from each data segment,
    while assuming that the slope in the right segment is actually 0.0.
    
    :Params x, y:
        The predictor and predictand datasets as lists of floats. The length of
        x and y must be the same size, but the values may contain an arbitrary
        number of missing values as specified by missing_val.
    :Param breakpoint:
        The index of the suspected changepoint in the data. Must be
        in range(len(x))
    :Param missing_val:
        A float or int representing the placeholder value for missing
        values in the dataset
    :Param vals:
        A 'Values' object for storing important intermediate values from the 
        computations during this test. This test requires vals to already 
        contain to values:
            vals.sse_resid - the sum square error on the residuals computed by
                fitting a linear regression on the entire dataset
            vals.n_vals - the number of valid (x, y) pairs given by the user
            vals.left_kth_regress - the regression on the left data segment
            vals.n_left - the number of valid (x, y) pairs in the left segment
            vals.right_kth_regress - the regression on the right data segment
            vals.n_right - the number of valid (x, y) pairs in the right segment
    :Returns:
        result, a dictionary with the following fields:
            :Param cmodel:
                The name of the test ('KTHTPR4')
            :Param bic: 
                The Bayes Information Criterion statistic for this analysis
            :Param mu: 
                A 2-value tuple with the y-intercepts of the regression here
            :Param alpha: 
                A 2-value tuple with the computed slopes of the regressions 
                here
            :Param sse_bic: 
                The sum square error on the computed BIC value
            :Param seg_lens: 
                A 2-value tuple with the length of the segments in the 
                regressions here
            :Param test_stat: 
                The test statistic value computed here
            :Param crit_val: 
                The critical value for the test statistic at alpha=0.05.
            :Param amp_est: 
                An estimate of the changepoint amplitude
        vals, the original 'Values' object provided to this method with no
            modifications.
                
    """
    cmodel = "KTHTPR4"
    
    left_x = x[:breakpoint+1]
    
    ## Recall the regression fit on the left and right segments
    kthl_left = vals.left_kth_regress
    kthl_right = vals.right_kth_regress
    
    y1 = kthl_left.y_int + kthl_left.slope*left_x[-1]
    y2 = kthl_right.y_med
    amp_est = y1 - y2
    
    count = kthl_left.nval + kthl_right.nval    
    
    # Compute test statistic, f1
    # Recall the residuals computed from KTHSLR1, where we assumed constant
    # slope throughout the data segment (sse_resid)
    sse_total = kthl_left.sseslope + kthl_right.sseflat
    f_val = ((vals.sse_resid-sse_total)/1.0)/(sse_total/(count-3))
    f_crit = lookup_critical(count, "f1")
    bic, bic_error, bic_dof_penalty = bayes(count, sse_total, 4)
    # output string
    head = "%7s %6.2f %7.2f %7.2f" % (cmodel, bic, bic_error, bic_dof_penalty)
    stats =  " %7.2f %7.2f %7.3f ------- %7.2f %7.2f" % (kthl_left.y_int, kthl_right.y_med,
                                                         kthl_left.slope, f_val, f_crit)
    tail = "% 7.2f %5d %4d" % (amp_est, vals.n_left, vals.n_right)
    #print (head+stats+tail)
    
    ## Record result of this analysis
    result = dict(cmodel=cmodel,
                  bic=bic, 
                  mu=[kthl_left.y_int, kthl_right.y_med],
                  alpha=[kthl_left.slope, 0.0],
                  sse_bic=sqrt(sse_total/vals.n_vals), # total number of non-missing values, from KTHSLR0
                  seg_lens=[vals.n_left, vals.n_right],
                  test_stat=f_val,
                  crit_val=f_crit,
                  amp_est=amp_est)
    
    return result