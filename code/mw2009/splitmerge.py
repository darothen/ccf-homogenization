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
from operator import itemgetter

# ccf-homogenization imports
from util import get_valid_data, compute_mean
from util import scale_series, compute_monthly_anomalies, imo2iym, within
from mw2009.chgptmodels import bayes, kth_line, t_test, lookup_critical

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

def splitmerge(network, beg_year=1, end_year=2, **kwargs):
    
    ## EXPERIMENTAL PLACEHOLDERS - will eventually be replaced with a master
    ## loop to do all the id pairs.
    
    #id1 = "215887"
    #id2 = "215615"
    id2 = "153430"
    id1 = "034572"
    
    ## 
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
                first_set = True
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
        print breakpoints
        print seg_bounds
        
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
            homog_segs = sorted(homog_segs, key=itemgetter(0))
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
        print new_breakpoints
        seg_bounds = zip(new_breakpoints[:-2], new_breakpoints[2:])
        print seg_bounds
        
        remove_breakpoints = set()
        merged_breakpoints = set()
        if iter > 0:
            
            print "Merge segments (isplit = 0), ipass: ", iter
                        
            for ((l ,r), new_bp) in zip(seg_bounds, new_breakpoints[1:-1]):
                
                print imo2iym(l), imo2iym(r)
                
                # short circuit and skip this segment if we already know that it's
                # homogenous
    #            this_within = lambda seg: within((l, r), seg)
    #            within_stable_segs = map(this_within, homog_segs)
    #            if any(within_stable_segs):
    #                print "Stable segment: ", imo2iym(l), imo2iym(r)
    #                if l == first: 
    #                    new_breakpoints.append(first)
    #                seg_lookup.append(((l, r), 'stable'))
    #                continue
        
        ## Apply the same adjustments and the same standard normal homogeneity
        ## test that we did in the previous splitting process. There is no 
        ## difference here until we consider what to do if we find a new 
        ## homogeneous segment.
                adjust = int(seg_bounds.index((l, r)) > 0)
                segment = diff_data[l+adjust:r+1]
                
                z = standardize(segment, MISS)
                likelihood_ratios = snht(z, MISS, standardized=True)
                z_count = len(get_valid_data(z))
                    
                ind_max_ratio = 0.0
                max_ratio = 0.0
                clip_ratios = likelihood_ratios[2:-2] # We clip the beginning and end
                for (ind, ratio) in zip(xrange(len(clip_ratios)), clip_ratios):
                    if ratio > max_ratio:
                        ind_max_ratio = ind
                        max_ratio = ratio
                        
                crit_val = lrt_lookup(z_count)
                possible_changepoint = l + ind_max_ratio + 2 + adjust
                
                y_new, m_new = imo2iym(possible_changepoint)
                
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
    print ym_breakpoints
    
    ## ENTERING MINBIC
    #left, bp, right = breakpoints[6:9]
    #left, bp, right = breakpoints[1:4]
    
    bp_dictionary = dict()
    for left,bp,right in zip(breakpoints[0:], breakpoints[1:], breakpoints[2:]):
        
        results_dictionary = dict()
        
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
        print "Entering MINBIC - %4d %2d    %4d %2d    %4d %2d" % (y1, m1, yb,
                                                                   mb, y2, m2)
        
        ## Print header for BIC changepoint testing -
        left_header = " QTYP     QVAL    QRSE     QPF     MU1     MU2  ALPHA1"
        right_header = "  ALPHA2   MSTAT   MCRIT    MOFF KNT1 KNT2"
        print (left_header+right_header)
        
        ################################################################################     
        
        ## Looks like the first test is KTHSLR1, kendall-theil method with slope calc
        ## We perform this test on the entire interval containing the breakpoint
        (seg_x, seg_data) = range(left_shift, right_shift+1), diff_data[left:right+1] 
        cmodel = "KTHSLR1"
        #lsql = least_squares(seg_x, seg_data, MISS)
        kthl = kth_line(seg_x, seg_data, MISS)
        
        nobs = right-left+1
        slpmed = kthl.slope
        yintmed = kthl.y_int
        sseredmed = kthl.sseslope
        nval = kthl.nval
        qoff = 0.0
        
        qslr1, rsq1, rsq2 = bayes(nobs, sseredmed, 2)
        # output string
        head = "%7s %6.2f %7.2f %7.2f" % (cmodel, qslr1, rsq1, rsq2)
        stats =  " %7.2f ------- %7.3f ------- ------- -------" % (yintmed, slpmed)
        tail = "% 7.2f %5d ----" % (qoff, nval)
        print (head+stats+tail)
        
        ## add to results dictionary
        results_dictionary[cmodel] = dict(q=qslr1, 
                                          rmu=[yintmed, yintmed],
                                          slpq=[slpmed, slpmed],
                                          sseq=sqrt(sseredmed/nval),
                                          mkfinq=[nval, 0],
                                          qstat=0.0,
                                          qcrit=99.0,
                                          qoff=qoff)
        
        ################################################################################     
        
        ## Now we begin the two-phase regressions, where we use the kendall-theil fit
        ## on the data at either side of the breakpoint.
        ##
        ## The first regression assumes that the segments have the same slope
        left_seg = range(left_shift, bp_shift+1)
        left_data = diff_data[left:bp+1]
        right_seg = range(bp_shift+1, right_shift+1)
        right_data = diff_data[bp+1:right+1]
        
        # kendall-theil method with 0 sloped segments
        cmodel = "KTHTPR0"
        
        kthl_left = kth_line(left_seg, left_data, MISS)
        kthl_right = kth_line(right_seg, right_data, MISS)
        
        left_y_med = kthl_left.y_med
        right_y_med = kthl_right.y_med
        q_off = left_y_med - right_y_med
        
        n_left = kthl_left.nval
        n_right = kthl_right.nval
        n_total = n_left + n_right
        
        stat_test = t_test(left_data, right_data, MISS)
        t_val = stat_test.t_val
        t_crit = lookup_critical(n_total-2, "t")
        
        sse_sum = kthl_left.sseflat + kthl_right.sseflat
        qslr1, rsq1, rsq2 = bayes(n_total, sse_sum, 3)
        # output string
        head = "%7s %6.2f %7.2f %7.2f" % (cmodel, qslr1, rsq1, rsq2)
        stats =  " %7.2f %7.2f ------- ------- %7.2f %7.2f" % (left_y_med, right_y_med,
                                                              t_val, t_crit)
        tail = "% 7.2f %5d %4d" % (q_off, n_left, n_right)
        print (head+stats+tail)
        
        ## add to results dictionary
        results_dictionary[cmodel] = dict(q=qslr1, 
                                          rmu=[left_y_med, right_y_med],
                                          slpq=[0.0, 0.0],
                                          sseq=sqrt(sse_sum/nval),
                                          mkfinq=[n_left, n_right],
                                          qstat=t_val,
                                          qcrit=t_crit,
                                          qoff=q_off)
        
        ################################################################################     
        
        ## The third regression tests for a step change with equal (constant) sloped
        ## segments
        cmodel = "KTHTPR1"
        
        all_data = diff_data[left:right+1]
        nobs = right+1-left
        left_data = diff_data[left:bp+1]
        right_data = diff_data[bp+1:right+1]
        
        ## 1) Compute the mean for *all* of the data
        all_valid_data = get_valid_data(all_data, MISS)
        all_mean = compute_mean(all_valid_data, valid=True)
        
        ## 2) use kendall-theil method with single slope
        # This method is slightly different than the kth_line() method above.
        # First, we get only the valid data, and we pair it with the natural ordering
        # of the data, i.e. 1, 2, 3...
        valid_all = all_valid_data
        n_all = len(valid_all)
        range_all = range(1, n_all+1)
        
        valid_left = get_valid_data(left_data, MISS)
        valid_right = get_valid_data(right_data, MISS)
        
        n_left, n_right = len(valid_left), len(valid_right)
        range_left = range(1, n_left+1)
        range_right = range(n_left+1, n_all+1)
        
        # Second, generate paired slopes for the first segment.
        nslp = 0
        r_temp = []
        for i in range(n_left-1):
            for j in range(i, n_left):
                if range_left[j] != range_left[i]:
                    nslp = nslp + 1
                    r_temp.append( (valid_left[j]-valid_left[i])/
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
                    r_temp.append( (valid_right[j]-valid_right[i])/
                                   (range_right[j]-range_right[i]) )
                    
        #Fourth, find the median slope from all the ones we computed
        islope = 1
        if not islope:
            r_slope = 0.0
        else:
            r_temp = sorted(r_temp)
            imed = (nslp - 1)/2
            if (nslp%2)==1: imed = imed+1 # offset by one to right if odd
            r_slope = r_temp[imed]
            
        print "slope, ic, imet: %7.2f %5d %5d" % (r_slope, nslp, imed)
        
        # Fifth, compute the first segment intercept, y-median - slope*x-median
        imed = (n_left - 1)/2
        if (n_left%2)==1: imed = imed+1
        range_med = range_left[imed]
        valid_left = sorted(valid_left)
        data_med = valid_left[imed]
        left_y_int = data_med-r_slope*range_med
        print "Seg1 - Xmed, Ymed, slope, Yint: %7.2f %7.2f %7.2f %7.3f" % (range_med, 
                                                                           data_med,
                                                                           r_slope,
                                                                           left_y_int)
        # BUG: Again in chgptmodel.kendalltheill(), there is a bug on line 2339. Starting
        #     here, we over-write the medians we found in both lists, and use the second
        #     segment for all our computations! I reproduce that behavior ehre by using the
        #     generic range_med and data_med values for rXmed and rYmed. Should we not
        #     care about different medians for each segment?
        # Sixth, compute the second segment intercept
        imed = (n_right - 2)/2
        if (n_right%2)==1: imed = imed + 1
        range_med = range_right[imed]
        valid_right = sorted(valid_right)
        data_med = valid_right[imed]
        right_y_int = data_med-r_slope*range_med
        print "Seg2 - Xmed, Ymed, slope, Yint: %7.2f %7.2f %7.2f %7.3f" % (range_med, 
                                                                           data_med,
                                                                           r_slope,
                                                                           right_y_int)
        
        # Seventh, we compute root mean square error of the residuals
        residuals = [MISS]*n_all # residuals of the fit
        fit = [MISS]*n_all       # fitted regression line
        valid_count = 0          # total number of non-missing values used
        r_sum_sqr_x = 0.0     
        r_sum_sqr_e = 0.0        # sum square of residuals
                                 # r_slope - slope of linear regression line
                                 # r_t - slope error
        for i in range(n_all):
            if all_data[i] != MISS:
                valid_count = valid_count + 1
                if valid_count < n_left:
                    y_int = left_y_int
                else:
                    y_int = right_y_int
                residuals[i] = (y_int + r_slope*(i+1)) - all_data[i]
                fit[i] = y_int + r_slope*(i+1)
                r_sum_sqr_e = r_sum_sqr_e + residuals[i]**2
                r_sum_sqr_x = r_sum_sqr_x + (float(i+1) - data_med)**2
        
        r_se_sqr = r_sum_sqr_e / (valid_count - 2)
        r_sb = sqrt(r_se_sqr / r_sum_sqr_x)
        r_t = r_slope / r_sb
        
        ############## END KENDALLTHIELL()
        r_mu = (left_y_int, right_y_int)
        r_alpha = r_slope
        SSE_red = r_sum_sqr_e
        
        # 3) Now it looks like we compute residuals for all our data
        r_residuals = [MISS]*nobs
        rssx = 0.0
        rsse = [0.0, 0.0]
        for k in range(nobs):
            if all_data[k] != MISS:
                if k < (bp+1-left):
                    ind = 0
                else:
                    ind = 1
                r_residuals[k] = all_data[k] - r_mu[ind] - r_alpha*(k+1)
                rsse[ind] = rsse[ind] + r_residuals[k]**2
                rssx = rssx + (float(k+1) - all_mean)**2
        
        # 4) We now have the squared error and are basically done!
        r_sum_sqr_tot = sum(rsse)
        
        ############## END KTHTPR1()
        
        # At this point, we have a few things - 
        #    r_mu - the y_intercepts of each segment
        #    r_alpha - the slope
        #    r_sum_sqr_tot - sum sqr total of residuals
        #
        # We now print out the info and compute critical values, BIC
        count = n_all
        left_count, right_count = n_left, n_right
        sseful = r_sum_sqr_tot # just computed
        # F-statistic
        f_val = ((sseredmed-sseful)/1.)/(sseful/(count-3))
        f_crit = lookup_critical(count-3, "f1")
        qslr1, rsq1, rsq2 = bayes(count, sseful, 4)
        # amplitude change estimate
        y1 = r_mu[0] + r_alpha * range_all[bp+1-left]
        y2 = r_mu[1] + r_alpha * range_all[bp+1-right]
        est = y1-y2
        # k, we have finished this god-awful changepoint
        # output string
        head = "%7s %6.2f %7.2f %7.2f" % (cmodel, qslr1, rsq1, rsq2)
        stats =  " %7.2f %7.2f %7.3f ------- %7.2f %7.2f" % (r_mu[0], r_mu[1],
                                                             r_alpha, f_val, f_crit)
        tail = "% 7.2f %5d %4d" % (est, n_left, n_right)
        print (head+stats+tail)
        
        ## add to results dictionary
        results_dictionary[cmodel] = dict(q=qslr1, 
                                          rmu=r_mu,
                                          slpq=[r_alpha, r_alpha],
                                          sseq=sqrt(sseful/nval),
                                          mkfinq=[n_left, n_right],
                                          qstat=f_val,
                                          qcrit=f_crit,
                                          qoff=est)
        
        ################################################################################     
        
        ## The fourth regression tests for a step change with any sloped segments, i.e.
        ## a full two phase regression. We borrow some of hte calculatiosn above, for
        ## simplicity' sake
        cmodel = "KTHTPR2"
        
        all_data = diff_data[left:right+1]
        nobs = right+1-left
        
        left_seg = range(left_shift, bp_shift+1)
        left_data = diff_data[left:bp+1]
        right_seg = range(bp_shift+1, right_shift+1)
        right_data = diff_data[bp+1:right+1]
        
        kthl_left = kth_line(left_seg, left_data, MISS)
        kthl_right = kth_line(right_seg, right_data, MISS)
        
        y1 = kthl_left.y_int + kthl_left.slope*left_seg[-1]
        y2 = kthl_right.y_int + kthl_right.slope*right_seg[0]
        est = y1 - y2
        
        count = len(get_valid_data(left_data, MISS))+len(get_valid_data(right_data, MISS))
        
        SSEful = kthl_left.sseslope + kthl_right.sseslope
        print kthl_left.sseslope, kthl_right.sseslope, sseredmed
        
        # note - sseredmed is the sse residuals of the slope from the kth_line for the 
        #    entire data segment
        f2_val = ((sseredmed-SSEful)/2.)/(SSEful/(count-4))
        f2_crit = lookup_critical(count-4, "f2")
        qslr1, rsq1, rsq2 = bayes(count, SSEful, 5)
        # output string
        head = "%7s %6.2f %7.2f %7.2f" % (cmodel, qslr1, rsq1, rsq2)
        stats =  " %7.2f %7.2f %7.3f %7.3f %7.2f %7.2f" % (kthl_left.y_int, kthl_right.y_int,
                                                           kthl_left.slope, kthl_right.slope,
                                                           f2_val, f2_crit)
        tail = "% 7.2f %5d %4d" % (est, n_left, n_right)
        print (head+stats+tail)
        
        ## add to results dictionary
        results_dictionary[cmodel] = dict(q=qslr1, 
                                          rmu=[kthl_left.y_int, kthl_right.y_int],
                                          slpq=[kthl_left.slope, kthl_right.slope],
                                          sseq=sqrt(SSEful/nval),
                                          mkfinq=[n_left, n_right],
                                          qstat=f2_val,
                                          qcrit=f2_crit,
                                          qoff=est)
        
        ################################################################################     
        
        ## The fifth regression tests for a step change with flat-to-sloped segments.
        ## Again, we will make use of some of the other values we have found.
        cmodel = "KTHTPR3"
        
        y1 = kthl_left.y_med
        y2 = kthl_right.y_int + kthl_right.slope*right_seg[0]
        est = y1 - y2
        
        count = len(get_valid_data(left_data, MISS))+len(get_valid_data(right_data, MISS))
        
        SSEful = kthl_left.sseflat + kthl_right.sseslope
        # note - sseredmed is the sse residuals of the slope from the kth_line for the 
        #    entire data segment
        f_val =  ((sseredmed-SSEful)/1.0)/(SSEful/(count-3))
        f_crit = lookup_critical(count, "f1")
        qslr1, rsq1, rsq2 = bayes(count, SSEful, 4)
        # output string
        head = "%7s %6.2f %7.2f %7.2f" % (cmodel, qslr1, rsq1, rsq2)
        stats =  " %7.2f %7.2f ------- %7.3f %7.2f %7.2f" % (kthl_left.y_med, kthl_right.y_int,
                                                            kthl_right.slope, f_val, f_crit)
        tail = "% 7.2f %5d %4d" % (est, n_left, n_right)
        print (head+stats+tail)
        
        ## add to results dictionary
        results_dictionary[cmodel] = dict(q=qslr1, 
                                          rmu=[kthl_left.y_med, kthl_right.y_int],
                                          slpq=[0.0, kthl_right.slope],
                                          sseq=sqrt(SSEful/nval),
                                          mkfinq=[n_left, n_right],
                                          qstat=f_val,
                                          qcrit=f_crit,
                                          qoff=est)
        
        ################################################################################     
        
        ## The sixth regression tests for a step change with sloped-to-flat segments.
        ## Again, we will make use of some of the other values we have found.
        cmodel = "KTHTPR4"
        
        y1 = kthl_left.y_int + kthl_left.slope*left_seg[-1]
        y2 = kthl_right.y_med
        est = y1 - y2
        
        count = len(get_valid_data(left_data, MISS))+len(get_valid_data(right_data, MISS))
        
        SSEful = kthl_left.sseslope + kthl_right.sseflat
        f_val = ((sseredmed-SSEful)/1.0)/(SSEful/(count-3))
        f_crit = lookup_critical(count, "f1")
        # note - sseredmed is the sse residuals of the slope from the kth_line for the 
        #    entire data segment
        qslr1, rsq1, rsq2 = bayes(count, SSEful, 4)
        # output string
        head = "%7s %6.2f %7.2f %7.2f" % (cmodel, qslr1, rsq1, rsq2)
        stats =  " %7.2f %7.2f %7.3f ------- %7.2f %7.2f" % (kthl_left.y_int, kthl_right.y_med,
                                                             kthl_left.slope, f_val, f_crit)
        tail = "% 7.2f %5d %4d" % (est, n_left, n_right)
        print (head+stats+tail)
        
        ## add to results dictionary
        results_dictionary[cmodel] = dict(q=qslr1, 
                                          rmu=[kthl_left.y_int, kthl_right.y_med],
                                          slpq=[kthl_left.slope, 0.0],
                                          sseq=sqrt(SSEful/nval),
                                          mkfinq=[n_left, n_right],
                                          qstat=f_val,
                                          qcrit=f_crit,
                                          qoff=est)
        
        ################################################################################ 
        ## We've finished the basic TPR tests. The PHA has loop-back criteria
        ## which could cause the code to double back and compute these test
        ## statistics over again or with subtle variations. I'll come back to
        ## that...
        ##
        ## For now, let's look at the test results and pick the most likely
        ## changepoint class
        
        results = [(cmodel, results_dictionary[cmodel]['q']) for cmodel in results_dictionary]
        results = sorted(results, key=itemgetter(1))
        cmodel, q = results[0]
        
        output = results_dictionary[cmodel]
        sseq=output['sseq']
        rmu=output['rmu']
        slpq=output['slpq']
        qstat=output['qstat']
        qcrit=output['qcrit']
        qoff=output['qoff']
        mkfinq=output['mkfinq']
        
        print ("Post: - - %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %5d %5d" %
               (q, sseq, rmu[0], rmu[1], slpq[0], slpq[1], qstat, qcrit, qoff,
                mkfinq[0], mkfinq[1]) )
        
        ## Store information about this breakpoint and its crucial stats
        bp_dictionary[bp] = dict(asigx=qoff, azscr=qoff/sseq,
                                 rslp=[kthl_left.slope, kthl_right.slope])
                              
        ## Also print an odd BIC: statement?
        # skip for now - basically prints the changepoint results, as well as
        # info about where the breakpoint is and the kthl slope on each side.

    ## Final stage - print the adjustment summaries
    sorted_bps = sorted(bp_dictionary.keys())
    for bp in sorted_bps:
        stats = bp_dictionary[bp]
        
        asigx=stats['asigx']
        azscr=stats['azscr']
        rslp=stats['rslp']
        
        end1 = bp
        y_end1, m_end1 = imo2iym(bp)
        beg2 = bp+1
        y_beg2, m_beg2 = imo2iym(bp+1)
        
        print ("%6s-%6s  --  -- MD TESTSEG ADJ: %7.2f %7.2f %8.4f %8.4f %5d %5d %3d %5d %5d %3d" % 
               (id1,id2, asigx, azscr, rslp[0], rslp[1], end1, y_end1, m_end1, beg2, y_beg2, m_beg2))
        
        
            
