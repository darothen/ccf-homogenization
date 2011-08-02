#!/usr/bin/env python
#
# Copyright (C) 2011 Daniel Rothenberg.
# See Google Code project page for license, 
# http://code.google.com/p/ccf-homogenization/

"""Port of splitmerge routine and utility functions 

"""
__docformat__ = "restructuredtext"

#http://docs.python.org/library/copy.html
from copy import deepcopy
# http://docs.python.org/library/math.html
from math import sqrt
# http://docs.python.org/library/operator.html
import operator
#http://docs.python.org/library/itertools.html
from itertools import combinations

import time

# ccf-homogenization imports
from util import get_valid_data, compute_mean, sign
from util import scale_series, compute_monthly_anomalies, imo2iym, within
from mw2009.chgptmodels import minbic

MINLEN = 18

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
    
    ## Find the valid data to use to compute the data_mean, etc.
    valid_data = get_valid_data(data, missing_val)
    num_vals = len(valid_data)
    data_mean = compute_mean(valid_data, valid=True)

    ## Compute the sum the squared error for each term (variance)
    variance_sum = 0.0
    for d in valid_data:
        variance_sum = variance_sum + (d-data_mean)**2

    ## The standard deviation is the root of the sum of the squared error
    sum_std = sqrt(variance_sum/(num_vals-2))
    
    ## Normalize each data value using this standard deviation
    standardized_data = []
    #for d in data[left:right+1]:
    for d in data:
        if d!= missing_val:
            standardized_data.append((d-data_mean)/sum_std)
        else:
            standardized_data.append(missing_val)
            
    return standardized_data

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
    
def snht(data, missing_val=-9999, valid_count=None, standardized=False):
    """Standard normal homogeneity test
    
    Loops over the given data and computes the likelihood ratio statistic for
    every value, using that value as the pivot to divide the data in two
    segments.
    
    :Param data:
        The data, as a list of floats, on which to conduct the test.
    :Param missing_val:
        (optional) The placeholder for missing values which should be excluded
        from computations.
    :Param valid_count:
        (optional) The number of valid values in the dataset.
    :Param standardized:
        (optional) Boolean flag indicating whether or not the data has already
        been standardized into a reference series. If not, will invoke the 
        standardization procedure on the data.
    :Return:
        A list with the same shape as the original list of data, containing
        the likelihood ratio test statistic for each data point. Missing values
        in the input data will be carried into the output list.
    
    """
    
    ## Standardize the data if this hasn't already been done
    if not standardized:
        data = standardize(data, missing_val)
        
    ## Return array of computed test statistics, initialized to missing_val's
    ts = [missing_val for d in data]
    
    if not valid_count:
        valid_count = len(get_valid_data(data, missing_val))

    pivot_count = range(valid_count-1)
    
    # BUG: This is *really* counter-intuitive, and probably a bug in the original 
    #     PHA code. Here, and in the PHA code, we use valid_count to effectively truncate
    #     the right tail of the data. The catch is, valid_count is the number of 'valid'
    #     data points - not *all* the data points. Because we end the right-seek
    #     at valid_count, we end up missing some of the data on the far right hand side
    #     of the array.
    #
    # -- CONFIRMED that this is a bug. Fix commented out below to be deployed
    #    when a new copy of the Fortran PHA is available.
    for pivot in pivot_count:
        
        ## Loop over the data, using each point as a pivot for computing the 
        ## likelihood ratio statistic.
        if data[pivot] != missing_val:
            left_series = get_valid_data(data[:pivot+1])
            right_series = get_valid_data(data[pivot+1:valid_count]) # BUG Line
            #right_series = get_valid_data(data[pivot+1:]) # FIX Line
            
            ## Compute mean of data left of the pivot; skip to next pivot value
            ## if no good data was found in this segment.
            sum_left = sum(left_series)
            nleft = len(left_series)
            if nleft != 0:
                mean_left = sum_left/nleft
            else:
                break
            
            ## Do the same for the data right of the pivot.
            sum_right = sum(right_series)
            nright = len(right_series)
            if nright != 0:
                mean_right = sum_right/nright
            else:
                break
            
            ## Compute and store test statistics
            #if pivot == valid_count-2:
            #    print nleft, mean_left, sum_left, nright, mean_right, sum_right
            ts[pivot] = nleft*(mean_left**2) + nright*(mean_right**2)
    
    return ts

def splitmerge(network, beg_year=1, end_year=2, **kwargs):
    
    ## EXPERIMENTAL PLACEHOLDERS - will eventually be replaced with a master
    ## loop to do all the id pairs.
    id_list = network.stations.keys()
    pair_results = dict()
    
    def dict_to_tuples(d):
        keys = d.keys()
        return [(key, d[key]) for key in keys]
    ## Generate station pairs for use in splitmerge by iteratively going through the
    ## station_list and adding stations in order of decreasing correlation. Skip a 
    ## neighbor if the pair is already present; want 20 stations or until all the
    ## correlated neighbors are used up.
#    pairs = []
#    for id1 in id_list:
#        neighbors = dict_to_tuples(network.correlations[id1])
#        sorted_neighbors = sorted(neighbors, key=operator.itemgetter(1))
#        added_pairs = 0
#        while sorted_neighbors and (added_pairs < 5):
#            id2, _ = sorted_neighbors.pop()
#            ordered_pair = tuple(sorted((id1, id2)))
#            if not ordered_pair in pairs:
#                pairs.append(ordered_pair)
#                added_pairs += 1
    import pickle
    pairs = list(pickle.load(open("fortran_pairs", "r")))
    #pairs = [('314684', '315177'),
    #         ('086997', '314684'), 
    #         ('307633', '314684'),
    #         ('124837', '324418'),
    #         ('164407', '314055')]
    #pairs = [('116738', '231037'), ]
    
    for (id1, id2) in pairs:
        print "Pair %s with %s" % (id1, id2)
        #continue
        pair_str = "%6s-%6s" % (id1, id2)
        #if pair_str != "830006-830009":
        #if pair_str != "215615-215887":
        #if pair_str != "111280-124837":
        #    continue
        
        raw_series = network.raw_series
        stations = network.stations
        series_copy = deepcopy(raw_series)
        
        min_ann = 5
        num_years = end_year - beg_year
        num_months = num_years*12
            
        for s in series_copy.itervalues():
            data = s.series
            scaled = scale_series(data, 0.1, s.MISSING_VAL)
            anomalies = compute_monthly_anomalies(scaled, s.MISSING_VAL)
            s.set_series(anomalies, s.years)
        
        ## Retrieve the data for each of the stations.
        station1 = stations[id1]
        series1 = series_copy[id1]
        data1 = series1.monthly_series
                
        station2 = stations[id2]
        series2 = series_copy[id2]
        data2 = series2.monthly_series
        
        
        #print data1[:50]
        #print data2[:50]
        #print "################################################################"
        ## Compute the difference series        
        diff_data = diff(data1, data2)
        MISS = series1.MISSING_VAL # Missing value placeholder
        
        ## Quickly pass through the data to find where it starts. We need to do this
        ## because it's possible that beg_year is earlier than the first year of 
        ## valid data in either data1 or data2. Furthermore, the original PHA code
        ## deliberately clipped off the first year of good data, so emulate that 
        ## effect here as well.
        ##
        ## Ultimately, we save the extreme early and extreme late month with valid
        ## data to use as our first guess at the undocumented changepoints.
        first = 0
        first_set = False
        last = 0
        for (i, d1, d2) in zip(xrange(num_months), data1, data2):
            if d1!=MISS and d2!=MISS:
                if first < 12:
                    first = i
                    #first_set = True
                #if not first_set:
                #    first = i
                #    first_set = True
                last = i
                
        ## Set the initial breakpoints and the list of already-found, homogenous
        ## segments.    
        breakpoints = [first, last, ]
        homog_segs = []
        
        #####################################################################
        ## BEGIN SPLITMERGE PROCESS TO GENERATE FIRST GUESS AT UNDOCUMENTED
        ## CHANGEPOINTS
        iter = 0 # counts how many times we've repeated the splitmerge process
        enter_BIC = False # break out of iterations into the BIC process?
        last_breakpoints = []
        while (iter < 10) and not enter_BIC:
            
            seg_bounds = zip(breakpoints[:-1], breakpoints[1:])
            last_breakpoints = deepcopy(breakpoints)
            new_breakpoints = deepcopy(breakpoints)
                
            new_homog_segs = []
        
            print "Parse segments (isplit = 1), ipass: ", iter
            #print breakpoints
            #print seg_bounds
            
            #####################################################################
            ## Loop over the segments we have just found and apply the 
            ## semi-hierarchical splitting algorithm and a 5% significance level
            ## to find where there is a split.
            for (l, r) in seg_bounds:
                
                y1, m1 = imo2iym(l)
                y2, m2 = imo2iym(r)
        
                # To be entirely faithful to the original Fortran PHA, we need to 
                # slightly adjust the left bound in this segment. Unless the left 
                # bound is the ultimate left bound on the difference series, then we
                # want to exclude it from the computations on this segment, since it
                # was included as a possible changepoint in the last segment. We
                # will use the variable adjust to force this slight modification.
                adjust = int(seg_bounds.index((l,r)) > 0)
                segment = diff_data[l+adjust:r+1]
                #numyr = r-l
                
            ## If there's less than min_ann years of data, then we can't handle the
            ## statistical analysis here and we will short-circuit the computations
            ## on this segment, proceeding to the next.
                if (r-l) <= min_ann:
                    print "Too short: ", imo2iym(l), imo2iym(r)
                    continue
                
            ## If we've previously found that this segment is homogenous (has no
            ## potential changepoint), then we can skip it as well and proceed to
            ## the next one.
                # Set the within() method to check if this segment is within any
                # previously found homogenous ones. Use lambda, since we can't pass
                # keyword or positional arguments to map().
                within_this_seg = lambda seg: within((l, r), seg)
                within_stable_segs = map(within_this_seg, homog_segs)
                if any(within_stable_segs):
                    print "Stable segment: ", imo2iym(l), imo2iym(r)
                    if l == first: 
                        new_breakpoints.append(first)
                    continue
                
            ## The standard normal homogeneity test - which is the statistical test
            ## we'll use to see if there is a potential changepoint in this segment
            ## - requires us to normalize our paired difference series. We can do
            ## that in snht(), but we'll do it right now so we can inspect those
            ## standardized values later.
                z = standardize(segment, MISS)
                #print segment[:50]
                #print z[:50]
                
            ## Apply standard normal homogeneity test. 
            ## For mechanics, see Alexandersson and Moberg 1997, Int'l Jrnl of
            ## Climatology (pp 25-34)
                likelihood_ratios = snht(z, MISS, standardized=True)
                z_count = len(get_valid_data(z))
                        
            ## We're left with the likelihood ratio for each value being a potential
            ## changepoint. Find the max ratio, and if that value is significant, let
            ## it be the newest potential changepoint.
                ind_max_ratio = 0
                max_ratio = 0.0
                clip_ratios = likelihood_ratios[2:-2] # clip the beginning and end,
                                                      # they can't be changepoints.
                for (ind, ratio) in zip(xrange(len(clip_ratios)), clip_ratios):
                    if ratio > max_ratio:
                        ind_max_ratio = ind
                        max_ratio = ratio
            ## Now we find the critical value for this data set, and check our max
            ## likelihood ratio against it
                crit_val = lrt_lookup(z_count)
                
                # The possible changepoint is the index of the max ratio we found. 
                # We have to shift it the following ways to align it to the original
                # data -
                #    1) shift by 2 re-aligns it from clip_ratios to likelihood_ratios
                #    2) shift by adjust re-aligns it to this segment in diff_data
                #    3) shift by l re-aligns it to the first index in diff_data
                possible_changepoint = l + ind_max_ratio + 2 + adjust
                
                y_new, m_new = imo2iym(possible_changepoint) # year, month
                
            ## If this is the first iteration, we indicate as such, and add the new
            ## changepoint
                if iter == 0: 
                    print "%6s-%6s MD        FIRST series %4d %2d to %4d %2d | at %4d %2d ts: %4.2f limit >: %3.2f" % (id1,id2,y1,m1,y2,m2,y_new,m_new,max_ratio,crit_val)
                    breakpoints.append(possible_changepoint)
                    breakpoints = sorted(breakpoints)
            
                else:
            ## Else, if we found a new possible changepoint, add it to our list.
                    if max_ratio > crit_val:
                        print "%6s-%6s MD Inhomogenity for series %4d %2d to %4d %2d | at %4d %2d ts: %4.2f limit >: %3.2f %4d" % (id1,id2,y1,m1,y2,m2,y_new,m_new,max_ratio,crit_val,z_count)
                        new_breakpoints.append(possible_changepoint)
                        
            ## If not, record that we found a homogeneous segment.   
                    else:
                        print "%6s-%6s MD      Homogeneous series %4d %2d to %4d %2d | at %4d %2d ts: %4.2f limit >: %3.2f %4d" % (id1,id2,y1,m1,y2,m2,y_new,m_new,max_ratio,crit_val,z_count)
                        new_homog_segs.append((l, r))
            
            ## Now we need to update our account of which segments were homogeneous,
            ## because we need to know during the next iteration. We will do this,
            ## as well as condense stable segments that lie adjacent to each other
            ## i.e, if we have the segments [(1,5), (5, 10,),, (12, 15)], then we 
            ## really have [(1,10), (12, 15)].
            homog_segs.extend(new_homog_segs)
            if homog_segs:
                homog_segs = sorted(homog_segs, key=operator.itemgetter(0))
                final_homog_segs = [homog_segs[0], ] # this will be like a stack
                for seg in homog_segs[1:]:
                    last_seg = final_homog_segs[-1]
                    if last_seg[1] == seg[0]:
                        new_seg = (last_seg[0], seg[1])
                        final_homog_segs.pop()
                        final_homog_segs.append(new_seg)
                    else:
                        final_homog_segs.append(seg)
                homog_segs = final_homog_segs
        
            ## So we have new segments that can be generated from these new
            ## breakpoints. Now, the PHA routine enters a "merge" process
            ## to see whether or not to keep these newly found changepoints or throw
            ## them out as false alarms. 
            ##
            ## We do this by "leapfrogging" every other breakpoint. This gives us
            ## a set of segments that all have another breakpoint in them. We want
            ## to see if these segments are homogeneous, because if they are, it
            ## means that the breakpoint we previously found in the segment has 
            ## been superseded.
            new_breakpoints = sorted(new_breakpoints)
            #print new_breakpoints
            seg_bounds = zip(new_breakpoints[:-2], new_breakpoints[2:])
            #print seg_bounds
            
            remove_breakpoints = set()
            merged_breakpoints = set()
            if iter > 0:
                
                print "Merge segments (isplit = 0), ipass: ", iter
                            
                for ((l ,r), new_bp) in zip(seg_bounds, new_breakpoints[1:-1]):
                    
                    #print imo2iym(l), imo2iym(r)
                    
    #                # short circuit and skip this segment if we already know that it's
    #                # homogenous
    #                this_within = lambda seg: within((l, r), seg)
    #                within_stable_segs = map(this_within, homog_segs)
    #                if any(within_stable_segs):
    #                    print "Stable segment: ", imo2iym(l), imo2iym(r)
    #                    if l == first: 
    #                        new_breakpoints.append(first)
    #                    seg_lookup.append(((l, r), 'stable'))
    #                    continue
                    # Set the within() method to check if this segment is within any
                    # previously found homogenous ones. Use lambda, since we can't pass
                    # keyword or positional arguments to map().
                    within_this_seg = lambda seg: within((l, r), seg)
                    within_stable_segs = map(within_this_seg, homog_segs)
                    if any(within_stable_segs):
                        print "Stable segment: ", imo2iym(l), imo2iym(r)
                        #if l == first: 
                        #    new_breakpoints.append(first)
                        merged_breakpoints.update([l, r])
                        continue
            
            ## Apply the same adjustments and the same standard normal homogeneity
            ## test that we did in the previous splitting process. There is no 
            ## difference here until we consider what to do if we find a new 
            ## homogeneous segment.
                    adjust = int(seg_bounds.index((l, r)) > 0)
                    segment = diff_data[l+adjust:r+1]
                    
                    z = standardize(segment, MISS)
                    likelihood_ratios = snht(z, MISS, standardized=True)
                    z_count = len(get_valid_data(z))
                        
                    ind_max_ratio = 0
                    max_ratio = 0.0
                    clip_ratios = likelihood_ratios[2:-2] # We clip the beginning and end
                    for (ind, ratio) in zip(xrange(len(clip_ratios)), clip_ratios):
                        if ratio > max_ratio:
                            ind_max_ratio = ind
                            max_ratio = ratio
                            
                    crit_val = lrt_lookup(z_count)
                    possible_changepoint = l + ind_max_ratio + 2 + adjust
                    
                    y_new, m_new = imo2iym(possible_changepoint)
                    
    
                    if z_count < 2:
                        y1, m1 = imo2iym(l)
                        y2, m2 = imo2iym(r)
                        print "%6s-%6s MD  No found peaks %4d %2d to %4d %2d" % (id1,id2,y1,m1,y2,m2)
                        print "%6s-%6s MD  Compress 1 out peak at %4d %2d" % (id1,id2,y_new,m_new)
                        #remove_breakpoints.add_
            ## If we found a new breakpoint that is statistically significant, then
            ## great! Let's keep it.
                    if max_ratio > crit_val:
                        print "%6s-%6s MD  Peak kept in merge at %4d %2d | ts: %4.2f limit >: %3.2f" % (id1,id2,y_new,m_new,max_ratio,crit_val)
                        merged_breakpoints.add(l)
                        merged_breakpoints.add(new_bp)
                        merged_breakpoints.add(r)
            ## If not, then this segment was homogeneous, so the breakpoint which
            ## already exists in it is no good.
                    else:
                        print "%6s-%6s MD Compress 2 out peak at %4d %2d | ts: %4.2f limit >: %3.2f" % (id1,id2,y_new,m_new,max_ratio,crit_val)
                        # Crap, if there are any potential breakpoints in this segment,
                        # we need to remove them because this segment is homogeneous. Let's
                        # remember this homogeneous segment for now and come back once
                        # we've found all of them.    
                        merged_breakpoints.update([l, r])
                        remove_breakpoints.add(new_bp)
            
            ## At this point, we have a set of all the breakpoints we've accumulated
            ## during this iteration of split/merge, as well as a set of breakpoints
            ## which we've found to be of no further use. We can difference update
            ## our set of breakpoints to remove these guys, and let those merged
            ## breakpoints be the set of newest breakpoints for the next splitmerge
            ## iteration.
                merged_breakpoints.difference_update(remove_breakpoints)
                breakpoints = list(merged_breakpoints)
            
            breakpoints = sorted(breakpoints)
            
            ## Did we actually find new breakpoints? If not, then we're done
            ## with splitmerge and can move on to the BIC process.
            enter_BIC = (breakpoints == last_breakpoints)
            iter = iter + 1
            
        ## Okay wow, we've potentially made it to the BIC stage now... !
        if first not in breakpoints:
            breakpoints.insert(0, first)
        ym_breakpoints = map(imo2iym, breakpoints)
        #print ym_breakpoints
        
        ## ENTERING MINBIC    
        bp_dictionary = dict()
####################################
##### MULTIPROCESS
        from multiprocessing import Pool
        

        global counter
        multi_bp_dict = {}
        counter = 0
        def cb(r):
            global counter
            #print counter, r
            counter += 1
        
        start = time.clock()         
        po = Pool(processes=4)
        for left,bp,right in zip(breakpoints[0:], breakpoints[1:], breakpoints[2:]):
                    
            if left != first:
                left = left + 1
            # recall that we only consider data after the first full year. we will be 
            # computing regressions with the independent variable indexed from this 
            # starting point, so we need to shift these indices. we also need to shift them
            # by +1 if this is any segment beyond the first one, so that we don't include
            # changepoints in more than one analysis.
            # TOTAL_SHIFT = -12 + 1 = -11
            # 
            # However, this shift is only necessary while looking at the array indices that
            # we generate using range(). the data should already be aligned correctly.
            total_shift = -12 + 1
            left_shift, bp_shift, right_shift = left+total_shift, bp+total_shift, right+total_shift
            y1, m1 = imo2iym(left)
            yb, mb = imo2iym(bp)
            y2, m2 = imo2iym(right)
            #print "Entering MINBIC - %4d %2d    %4d %2d    %4d %2d" % (y1, m1, yb,
            #                                                           mb, y2, m2)
            (seg_x, seg_data) = range(left_shift, right_shift+1), diff_data[left:right+1]
            bp_index = bp-left
            #print len(seg_x), len(seg_data), bp_index
            #bp_analysis = minbic(seg_x, seg_data, bp_index, MISS)
            multi_bp_dict[bp] = po.apply_async(minbic,(seg_x,seg_data,bp_index,MISS,),callback=cb)
        po.close()
        po.join()
        for bp in multi_bp_dict:
            r = multi_bp_dict[bp]
            multi_bp_dict[bp] = r.get()
        #print "counter - %d" % counter
        elapsed = (time.clock() - start)
        print "ELAPSED TIME - %2.3e" % elapsed
        #print new_bp_dict
####################################
##### NORMAL        
#        start = time.clock()
#        for left,bp,right in zip(breakpoints[0:], breakpoints[1:], breakpoints[2:]):
#                    
#            if left != first:
#                left = left + 1
#            # recall that we only consider data after the first full year. we will be 
#            # computing regressions with the independent variable indexed from this 
#            # starting point, so we need to shift these indices. we also need to shift them
#            # by +1 if this is any segment beyond the first one, so that we don't include
#            # changepoints in more than one analysis.
#            # TOTAL_SHIFT = -12 + 1 = -11
#            # 
#            # However, this shift is only necessary while looking at the array indices that
#            # we generate using range(). the data should already be aligned correctly.
#            total_shift = -12 + 1
#            left_shift, bp_shift, right_shift = left+total_shift, bp+total_shift, right+total_shift
#            y1, m1 = imo2iym(left)
#            yb, mb = imo2iym(bp)
#            y2, m2 = imo2iym(right)
#            print "Entering MINBIC - %4d %2d    %4d %2d    %4d %2d" % (y1, m1, yb,
#                                                                       mb, y2, m2)
#            (seg_x, seg_data) = range(left_shift, right_shift+1), diff_data[left:right+1]
#            bp_index = bp-left
#            #print len(seg_x), len(seg_data), bp_index
#            bp_analysis = minbic(seg_x, seg_data, bp_index, MISS)
#            
#            bp_dictionary[bp] = bp_analysis    
#        elapsed2 = (time.clock() - start)
#        print "ELAPSED TIME = %3.2e" % elapsed2
        
        ##################################3
        ## Print the adjustment summaries
        bp_dictionary = multi_bp_dict
        sorted_bps = sorted(bp_dictionary.keys())
        ndelete = []
        valid_bps = {}
        for bp in sorted_bps:
            stats = bp_dictionary[bp]
            
            cmodel=stats['cmodel']
            iqtype=stats['iqtype']
            asigx=stats['offset']
            azscr=stats['offset_z']
            rslp=stats['slopes']
            
            end1 = bp
            y_end1, m_end1 = imo2iym(end1)
            beg2 = bp+1
            y_beg2, m_beg2 = imo2iym(beg2)
            
            # If cmodel is *SLR*, then there is no breakpoint
            if 'SLR' in cmodel:
                print ("%s-%s  --  -- MD TESTSEG SKIP: %7.2f %5d %5d %3d %5d %5d %3d" %
                       (id1, id2, asigx, end1, y_end1, m_end1, beg2, y_beg2, m_beg2))
                # Don't store it!
            else:
                print ("%6s-%6s  --  -- MD TESTSEG ADJ: %7.2f %7.2f %8.4f %8.4f %5d %5d %3d %5d %5d %3d %2d" % 
                       (id1,id2, asigx, azscr, rslp[0], rslp[1], end1, y_end1, m_end1, beg2, y_beg2, m_beg2, iqtype))
                # Store it!
                valid_bps[bp] = stats
        
        ###############################
        ## Go back and see if we can get rid of some of the change points.
        ## If 2 or more of the chgpts are within MINLEN,
        ##    a) if the chgpt estimates are the same sign, then test each
        ##        singly with same endpoints and keep lowest BIC 
        ##    b) if not the same sign,
        ##        retain earliest changepoint
        # add the first, last to valid_bps
        interior_bps = valid_bps.keys()
        # Add first, last if not already in interior_bps
        for bp in [first, last]:
            if bp not in interior_bps:
                interior_bps.append(bp)
        sorted_bps = sorted(interior_bps)
        for left in sorted_bps:
            print sorted_bps, left
            ## We're looking for the next interim breakpoint that satisfies two
            ## conditions:
            ##    1) at least MINLEN valid data (non-missing to the right)
            ##    2) has at least one breakpoint between 'left' and it
            right = 0
            close_bps = []
            for right in sorted_bps: 
                if right <= left: continue
                
                if not close_bps:
                    close_bps.append(right)
                else:
                    valid_between_bps = diff_data[close_bps[-1]:right]
                    valid_length = len(get_valid_data(valid_between_bps, MISS))
                    print imo2iym(close_bps[-1]),valid_length,imo2iym(right)
                    if valid_length > MINLEN:
                        break
                    close_bps.append(right)
            # We could actually run out of things in sorted_bps, and wind up with
            # right == close_bps[-1]. Detect that and break out of this analysis
            # if that happens.
            if close_bps[-1]==right: break
            
            if left != first:
                left = left + 1
            close_bp_results = {}
            for bp in close_bps:
                        
#                # recall that we only consider data after the first full year. we will be 
#                # computing regressions with the independent variable indexed from this 
#                # starting point, so we need to shift these indices. we also need to shift them
#                # by +1 if this is any segment beyond the first one, so that we don't include
#                # changepoints in more than one analysis.
#                # TOTAL_SHIFT = -12 + 1 = -11
#                # 
#                # However, this shift is only necessary while looking at the array indices that
#                # we generate using range(). the data should already be aligned correctly.
                total_shift = -12 + 1
                left_shift, bp_shift, right_shift = left+total_shift, bp+total_shift, right+total_shift
                y1, m1 = imo2iym(left)
                yb, mb = imo2iym(bp)
                y2, m2 = imo2iym(right)
                
                print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
                print y1,m1,"-",yb,mb,"-",y2,m2
                print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
                
                (seg_x, seg_data) = range(left_shift, right_shift+1), diff_data[left:right+1]
                bp_index = bp-left
                bp_analysis = minbic(seg_x, seg_data, bp_index, MISS, kthslr0_on=True)
                
                cmodel=bp_analysis['cmodel']
                iqtype= bp_analysis['iqtype']
                offset= bp_analysis['offset']
                rslp= bp_analysis['slopes']
                crit_val = bp_analysis['crit_val']
                test_stat = bp_analysis['test_stat']
                bic = bp_analysis['bic']
                
                print ("Interim chgpt: %s %4d %2d %4d %2d %4d %2d %8.2f %8.2f %8.2f %8.2f %7.3f %7.3f %2d" %
                       (pair_str, y1, m1, yb, mb, y2, m2, bic, test_stat, crit_val, offset, rslp[0], rslp[1], iqtype))                    
                
                close_bp_results[bp] = bp_analysis

            # Now we have a small problem... we might have more than one breakpoint,
            # so we need to choose which one is best. We will check the sign of
            # the breakpoint amplitude changes:
            sign_of_amps = map(sign, [close_bp_results[bp]['offset'] for bp in close_bps])
            positive = lambda x: sign(x) >= 0
            negative = lambda x: sign(x) <= 0
            zero = lambda x: sign(x) == 0
            print "------------>",[close_bp_results[bp]['offset'] for bp in close_bps]
            if (all(map(positive, sign_of_amps)) or 
                all(map(negative, sign_of_amps))):    
                # Pick the best (minimum BIC)          
                bics = [(bp, close_bp_results[bp]['bic']) for bp in close_bps]
                sorted_bics = sorted(bics, key=operator.itemgetter(1))
                smallest_bp = sorted_bics[0][0]
                
                # Remove this smallest-bic bp from the in-interval bps 
                close_bps.remove(smallest_bp)
                valid_bps[smallest_bp] = close_bp_results[smallest_bp] 
                
                #print "leftovers",close_bps
                for bp in close_bps: # The remaining bps which we will reject
                    sorted_bps.remove(bp) # Remove them from this loop
                    del valid_bps[bp] # Remove them as valid 
                    
                yb, mb = imo2iym(smallest_bp)
                print ("Same domain - Lowest Interim: %s %4d %2d" % 
                       (pair_str, yb, mb))
            elif (all(map(zero, sign_of_amps))):
                # Choose the earliest changepoint; the rest of these have
                # amplitude changes which are 0.
                first_bp, last_bp = close_bps[0], close_bps[-1]
                
                # Remove the first interim bp and update valid_bps with this new
                # computation. 
                close_bps.remove(first_bp)
                valid_bps[first_bp] = close_bp_results[first_bp]
                
                # Reject remaining interim bps
                for bp in close_bps:
                    sorted_bps.remove(bp)
                #    del valid_bps[bp]
                    
                yb, mb = imo2iym(first_bp)
                print ("Null domain - Earliest Interim : %s %4d %2d" %
                       (pair_str, yb, mb))
            else:
                # We'll use the earliest interim changepoint, but we need
                # to get rid of bad data. Replace all the data between the 
                # interim changepoints as missing and re-compute BIC.
                first_bp, last_bp = close_bps[0], close_bps[-1]
                first_bp_index = first_bp-left
                last_bp_index = last_bp-left
                
                print len(seg_x), len(seg_data)
                print first_bp_index+1, last_bp_index+1
                print left, bp, right
                for i in range(first_bp_index+1, last_bp_index+1):
                    print i
                    seg_x[i] = MISS
                    seg_data[i] = MISS
                    ndelete.append(i)
                bp_analysis = minbic(seg_x, seg_data, first_bp_index, MISS, kthslr0_on=True)
                
                # Remove the first interim bp and update valid_bps with this new
                # computation. 
                close_bps.remove(first_bp)
                valid_bps[first_bp] = bp_analysis
                
                # Reject remaining interim bps
                for bp in close_bps:
                    sorted_bps.remove(bp)
                #    del valid_bps[bp]
                
                yb, mb = imo2iym(first_bp)
                print ("Diff domain - Earliest Interim : %s %4d %2d" %
                       (pair_str, yb, mb))                
    
        valid_bps['del'] = ndelete
        pair_results[pair_str] = valid_bps
        
        #print "ELAPSED TIMES = %3.2e %3.2e" % (elapsed1, elapsed2)
    print "done"
    ##
    import pickle
    f = open("pair_results", 'w')
    pickle.dump(pair_results, f)
    return pair_results
            
            
