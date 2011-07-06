#!/usr/bin/env python
#
# Copyright (C) 2011 Daniel Rothenberg.
# See Google Code project page for license, 
# http://code.google.com/p/ccf-homogenization/

"""Port of utility functions used by splitmerge.v21.f 

"""
__docformat__ = "restructuredtext"

# http://docs.python.org/library/math.html
from math import sqrt
# ccf-homogenization imports
from util import get_valid_data, compute_mean

def diff(data1, data2, missing_val=-9999):
    """Computes the difference series between data1 and data2.
    
    diff[i] = data1[i] - data2[i] if neither data1[i] or data2[i] are
    equal to missing_val. If either one is, then diff[i] = missing_val.
    
    :Params data1, data2:
        The two lists of data which will be differenced against each
        other. Note that data1 and data2 must have the same length!
    :Param missing_val:
        The placeholder for missing values in the dataset.
    :Return:
        A list of data containing the difference series between data1
        and data2, with the same length as data1 and data2.
    
    """
    assert len(data1) == len(data2)
    
    diff_data = []
    for (d1, d2) in zip(data1, data2):
        if d1 != missing_val and d2 != missing_val:
            append_val = (d1-d2)
        else:
            append_val = missing_val
        diff_data.append(append_val)
        
    return diff_data

def standardize(data, missing_val=-9999):
    """Standardize a subset of data using a normal score.
    
    Computes the normal score for a set of data necessary for performing
    the standard normal homogenity test. The normal scores are defined as
    
    Z_i = (Q_i - Qbar) / sigma_Q,
    
    For more information, please refer to Alexandersson and Moberg, 1997,
    Int'l Journal of Climatology Vol. 17, pp 25-34.    
    
    This is a direct port of splitmerge.v21f.f > subroutine 'standard'.
    
    :Param data:
        The dataset to standardize.
    :Param missing_val:
        The placeholder for missing data.
    :Return:
        A list of length (right-left), with the standardized reference 
        values computed here.
    
    """
    
    ## Find the valid data to use to compute the mean, etc.
    #valid_data = get_valid_data(data[left:right+1], missing_val)
    valid_data = get_valid_data(data, missing_val)
    rNum = len(valid_data)
    rMean = compute_mean(valid_data, valid=True)
    rSum = sum(valid_data)

    ## Compute the sum the squared error for each term
    rVarSumm = 0.0
    for d in valid_data:
        rVarSumm = rVarSumm + (d-rMean)**2

    ## The standard deviation is the root of the TSE
    rsqVar = sqrt(rVarSumm/(rNum-2))
    
    ## Normalize each data value using this standard deviation
    rStd = []
    #for d in data[left:right+1]:
    for d in data:
        if d!= missing_val:
            rStd.append((d-rMean)/rsqVar)
        else:
            rStd.append(missing_val)
            
    #print rMean, rSum, rNum, rVarSumm, rsqVar
            
    return rStd

def lrt_lookup(num_vals):
    """Looks up the test statistic critical value corresponding to a significance
    level of 95% for the likelihood ratio test as formulated by Alexandersson
    and Moberg, given the number of data values used in the test.
    
    In Alexandersson and Moberg 1997, Int'l Jrnl of Climatology (pp 25-34),
    a test based on the likelihood ratio test is given to help to detect
    inhomogeneities in timeseries datasets. This function uses their lookup 
    table of critical values for the test from Appendix 2, Table AI of their
    paper, and linearly interpolated between these values to determine 
    non-specified critical values. 
    
    :Param num_vals:
        The number of values used in the statistical test.
    :Return:
        The critical value of the test statistic corresponding to that number
        of values.
        
    """
    sig95 = [4.54,5.70,6.95,7.65,8.10,8.45,8.65,8.80,8.95,9.05,
             9.15,9.35,9.55,9.70]
    nvals = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250]
    
    if num_vals < nvals[0]: 
        return 99999.0 # Too few values!
    elif num_vals >= nvals[-1]:
        return sig95[-1] # This is a good approximation as n->infinity
    else:
        # Select the two values from nsig which bound num_vals
        bounds = zip(nvals[:-1], nvals[1:])
        bi = 0
        
        l, r = bounds[bi]
        while not (l <= num_vals <= r): 
            bi = bi+1
            l, r = bounds[bi]
        left_ind, right_ind = nvals.index(l), nvals.index(r)
        
        # Estimate the critical value using linear interpolation between
        # the two lookup values we found
        bound_left, bound_right = l, r
        crit_left, crit_right = sig95[left_ind], sig95[right_ind]
        
        interp = num_vals - bound_left
        slope = (crit_right - crit_left) / (bound_right - bound_left)
        
        crit_val = crit_left + slope*interp
        return crit_val
    
def snht(data, missing_val=-9999, mcnt=None, standardized=False):
    """Standard normal homogeneity test
    
    """
    
    if not standardized:
        data = standardize(data, missing_val)
        
    # return array of computed test statistics, initialized to missing_val
    ts = [missing_val for d in data]
    
    if not mcnt:
        mcnt = len(get_valid_data(data, missing_val))

    pivot_count = range(mcnt-1)
    
    # BUG: This is *really* counter-intuitive, and probably a bug in the original 
    #     PHA code. Here, and in the PHA code, we use mcnt to effectively truncate
    #     the right tail of the data. The catch is, mcnt is hte number of 'valid'
    #     data points - not *all* the data points. Because we end the right-seek
    #     at mcnt, we end up missing some of the data on the far right hand side
    #     of the array.
    for pivot in pivot_count:
        
        if data[pivot] != missing_val:
            left_series = get_valid_data(data[:pivot+1])
            right_series = get_valid_data(data[pivot+1:mcnt])
            
            sum_left = sum(left_series)
            nleft = len(left_series)
            if nleft != 0:
                mean_left = sum_left/nleft
            else:
                break
        
            sum_right = sum(right_series)
            nright = len(right_series)
            if nright != 0:
                mean_right = sum_right/nright
            else:
                break
            
            ts[pivot] = nleft*(mean_left**2) + nright*(mean_right**2)
    
    return ts
                