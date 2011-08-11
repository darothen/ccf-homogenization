#!/usr/bin/env python
#
# Copyright (C) 2011 Daniel Rothenberg.
# See Google Code project page for license, 
# http://code.google.com/p/ccf-homogenization/
"""Final stage of MW2009 process - Estimation of adjusted series

"""
__docformat__ = "restructuredtext"

from util import imo2iym, get_valid_data, compute_monthly_anomalies
from util import compute_first_diff, compute_corr
from splitmerge import diff

from mw2009.chgptmodels import minbic, kthtpr0

import numpy as np
import operator

def estamt(network, minlenshf=24, **hom_params):
    """
    COPIED FROM ucpmonthly.v24a.f:
    
    The major steps in determining the best adjustment value for each station
    and changepoint. Entire network undergoes each of the following processes.
    In order:
        1) Remove unusable data. Align move swith respect to non-missing data
        and compress out changes that are too close AND the data between them.
        
        2) ISTEP=2 processing begins the adjustment process by removing the 
        non-significant changepoints to lengthen segments.
        
        3) NPASS (:= ISTEP=3) finishes the adjustment process by testing for the
        minimum number of months in a segment and number of neighbors with which
        the difference series can be examined.
        
        4) Final adjusted output is written.
    """
    
    ## FILTER 4
    ## Since the amplitude estimate MUST rely upon a minimum of MINLEN months to
    ## get even close to a reliable estimate at this point, it is assumed that
    ## the changepoints are as good as the station history files. Therefore,
    ## align moves with respect to non-missing data and compress out changes
    ## that are too close AND the data between them (i.e., less than MINLEN
    ## apart)
    
    #station_list = network.stations.keys()
    all_station_list = network.stations.keys()
    #station_list = ["215887", ]
    station_list = all_station_list
    
    # for each station...
    for id in station_list:
        station_index = station_list.index(id)
        station_series = network.raw_series[id]
        station_data = station_series.monthly_series[:]
        missing_val = station_series.MISSING_VAL
            
        # ... gen arrays for alignment
        move, amt, mday = [], [], []
        changepoints = station_series.changepoints
        cps = sorted(changepoints.keys())
        for cp in cps: 
            print "  Hist move: ",len(move)+1,station_index+1,imo2iym(cp) 
            move.append(cp)
            amt.append(changepoints[cp]['jsum'])
            mday.append(31)
        movnum = len(changepoints)
    
        if movnum > 0:
            ## At this point, the Fortran code executes alignmoves() in
            ## SHAPinp.v6c.f to reconcile the fact that station history files
            ## report dates of moves. It also removes segments that eare too short
            ## - less than minlenshf. Instead of implementing alignmoves(). Right
            ## now, I'll only implement this second functionality.
            #alignmoves()
            
            ####################################################################            
            # Seek to find first and last month indices
            first_set = False
            for month in range(len(station_data)):
                # Skip first year
                if month < 12: continue
                
                if station_data[month] != missing_val:
                    if not first_set:
                        first = month
                        first_set = True
                    last = month
                                
            cps = sorted(changepoints.keys())
            cps.insert(0, first)
            cps.append(last)
            for (cp1, cp2) in zip(cps[:], cps[1:]):
                if (cp2 - cp1) < minlenshf:
                    months_to_delete = range(cp1+1, cp2+1)
                    network.raw_series[id].delete_months(months_to_delete)
                    
                    if cp2 == last:
                        del_key = cp1
                    else:
                        del_key = cp2
                    if del_key in network.raw_series[id].changepoints:
                        #print len(network.raw_series[id].changepoints.keys()),
                        del network.raw_series[id].changepoints[del_key]
                        #print len(network.raw_series[id].changepoints.keys()),
                        #raw_input("pause")
                    
                    del_str = "Del 1st segment: " if cp1 == first else "Delete segment: "
                    
                    print id,station_index+1,del_str,imo2iym(cp1),cp1,imo2iym(cp2),cp2            
            
            new_changepoints = network.raw_series[id].changepoints
            new_cps = sorted(new_changepoints.keys())
            print "  First data value: ",imo2iym(first)
            for cp in new_cps:
                print "  End seg:",new_cps.index(cp)," ym: ",imo2iym(cp),cp,new_changepoints[cp]['jsum']
            print "    End segment ym: ",imo2iym(last),last
            
            # Finally, add first and last value to the list of changepoints.
            first_stats = dict(ahigh=0.0, astd=0.0, jsum=0)
            last_stats = dict(ahigh=0.0, astd=0.0, jsum=0)
            network.raw_series[id].changepoints[first] = first_stats
            network.raw_series[id].changepoints[last] = last_stats
            ####################################################################
    
    ## Series of debug print statements summarizing the final list of 
    ## changepoints. Not necessary at the moment 
                
    ############################################################################
    # The subnetwork processing became a multi-step process plus a "post-process
    # pass" to manage:
    #    1) problems with documented changepoints with NO undocumented support
    #    2) determine the best amplitude estimation for each confirmed changepoint
    
    for step in [2, 3]:
        ## Setup output strings based on the step used. Only cosmetic differences
        ## really.
        
        iminlen = hom_params['minlen']
        numclim = 3
        ## STEP 1 - NEVER USED (technically the history-consideration done previously
        if step == 1:
            continue
        elif step == 2:        
            ## STEP 2 - NOT SIG REMOVAL
            ## equivalent to ipass loopback for istep == 2 in Fortran PHA
            print " ---------------- NOT SIG REMOVAL --------------- "
            tstr = "Not sig: "
            outid = "NS"
            ipass = 1
        elif step == 3:
            ## STEP 3 - ADJUSTMENT OF DISCONTINUITIES
            # equivalent to ipass loopback for istep == in FORTRAN PHA
            print " ---------------- ADJUST DISCONTINUITY STEP --------------- "
            print "Adjpass, iminlen, numclim","--",iminlen,numclim
            print " ---------------- NPASS --------------- "
            tstr = "Dstep Dtrend: "
            outid = "WM"
            ipass = ipass + 1
        

        final_results = dict()
        print "  NET   STN    FILT TECH      ------ AFTER ------    ------ BEFORE ------"
        # Process each station and its network of neighbors
        for id in station_list:
            station_index = station_list.index(id)
            
            station_cp_dict = network.raw_series[id].changepoints 
            sorted_cps = sorted(station_cp_dict.keys())
                
            station_series = network.raw_series[id]
            missing_val = station_series.MISSING_VAL
            
            # compute monthly anomalies for this station data
            station_anomalies = station_series.monthly_anomaly_series
            
            # What are the first and last valid months in this station's data set?
            # We've saved them as the first and last changepoint before...
            first = sorted_cps[0]
            last = sorted_cps[-1]
    
            # What are the pairs to this station that we need to consider? 
            station_pairs = []
            for other_id in all_station_list:
                pair = tuple(sorted([id, other_id]))
                if pair in hom_params['pairs']:
                    station_pairs.append(pair)
            print station_pairs
    
            # List the remaining changepoints after the "confirmfilt" process
            for cp in sorted_cps:
                cp_stats = station_cp_dict[cp]
                hit_count = cp_stats['jsum']            
                iy, im = imo2iym(cp)
                print ("%3d %5d %6s Estamt chgin: -- %4d %2d %4d %3d" %
                       (ipass, station_index, id, iy, im, cp, hit_count) )
                
            ## ACCUMULATE PAIRED CHANGEPOINTS AND AMPLITUDE ESTIMATES
            # Loop over "brackets" of changepoints - that is, for changepoints
            # [a, b, c, d], consider the two brackets [a,b,c] and [b,c,d] with 
            # the center value of the changepoints. Note that in the Fortran PHA,
            # we go through these brackets in reverse order - right to left.
            brackets = zip(sorted_cps[-3::-1], sorted_cps[-2::-1], sorted_cps[::-1])
            final_results[id] = dict()
            for bracket in brackets:
            #for bracket in brackets[:1]:
                (left, cp, right) = bracket[:]
                
                ly, lm = imo2iym(left)
                cpy, cpm = imo2iym(cp)
                ry, rm = imo2iym(right)
                print "Oriented: ",'--','--','--',left,cp,cp+1,right
                
                # setup the output string for this bracket's tests
                chgptstr = ("  Win1: %5d %4d%2d %5d %4d%2dto Win2: %5d %4d%2d %5d %4d%2d" %
                            (left,ly,lm,cp,cpy,cpm,cp,cpy,cpm,right,ry,rm))
                
                ## THIS SECTION ACCUMULATES TARGET-NEIGHBOR COMPARISONS
                # See if there are enough homogeneous data in the target;
                # check each window
                valid_count_right = len(get_valid_data(station_data[cp+1:right+1], 
                                                       missing_val))
                valid_count_left = len(get_valid_data(station_data[left:cp+1],
                                                      missing_val))
                # if the segment length (valid count) is too short, skip this
                # changepoint (for now)
                if (valid_count_left < iminlen):
                    print "Adjpass seg2 short ",station_index,id,chgptstr,valid_count_left
                    continue
                if (valid_count_right < iminlen):
                    print "Adjpass seg1 short ",station_index,id,chgptstr,valid_count_right
                    continue
                
                ## We've pass the too-little-data pitfall. Now, we are actually going
                ## to go back through our paired neighbors and compute some final
                ## statistics about these changepoints. We'll store them in a
                ## dictionary for later, just like the pair_results dictionary
                ## from splitmerge
                pair_results = dict()
                #for (id1, id2) in [("215887", "200779")]:
                for (id1, id2) in station_pairs:
                    # Reset the left, cp, and right indices to the original
                    # bracket we're considering. We are going to be changing them
                    # while we look at this pair
                    (left, cp, right) = bracket[:]
                    
                    ## Figure out which station is the neighbor (not the target
                    ## we're currently considering). At the same time, note that if
                    ## the target is the 2nd changepoint, the adjustments will be
                    ## flipped in sign, so we need to have a correction factor ready
                    correction = 1.0
                    if id == id1:
                        neighb_id = id2
                    else:
                        neighb_id = id1
                        #correction = -1.0                    
                    
                    # Add this pair to pair_results if it's not already there
                    (ida, idb) = sorted([id1,id2])
                    pair_str = "%s-%s" % (ida, idb)
                    if pair_str not in pair_results:
                        pair_results[neighb_id] = dict()
                    print pair_str
                    
                    neighb_index = all_station_list.index(neighb_id)
                    neighb_cp_dict = network.raw_series[neighb_id].changepoints 
                           
                    neighb_series = network.raw_series[neighb_id]
                    neighb_anomalies = neighb_series.monthly_anomaly_series
    
                    ## Generature a difference data set for this pair of stations
                    diff_data = diff(station_anomalies, neighb_anomalies)
                    
                    ## It's possible that in the [left, right] bracket we're looking
                    ## at, there's a changepoint in the paired neighbor. We need
                    ## to adjust the endpoints of the bracket to exclude those
                    ## breakpoints
                    # Check right-hand side first and break out if ...
                    right_seg_len = len(get_valid_data(diff_data[cp+1:right+1]))
                    #right_seg_len = len(diff_data[cp+1:right+1])
                    for month in range(cp+1, right+1):
                        if month == last: continue
                        
                        # ... we hit a changepoint in the neighbor ...
                        if month in neighb_cp_dict:
                            neighb_hits = neighb_cp_dict[month]['jsum']
                            right_seg_len = len(get_valid_data(diff_data[cp+1:month+1]))
                            #right_seg_len = len(diff_data[cp+1:month+1])
                            print ("CHG2: ",neighb_index,neighb_id,"num,edit,2b,2e,imo,nhits",
                                   right_seg_len,"--",cp+1,right,month,neighb_hits)
                            
                            right = month
                            break
                    # ... and the final right-segment is too short
                    print left, cp, right
                    if right_seg_len < iminlen:
                        print ("Low2: ",neighb_index,neighb_id,"num,edit,2b,2e,imo,nhits",
                                right_seg_len,"--",cp+1,right,month,"--")
                        continue
                    
                    # Now, check the left-hand side and break out if ...
                    left_seg_len = len(get_valid_data(diff_data[left:cp+1]))
                    for month in range(cp-1, left, -1):
                        if month == first: continue
                        
                        # ... we hit a changepoint in the neighbor ...
                        if month in neighb_cp_dict:
                            neighb_hits = neighb_cp_dict[month]['jsum']
                            left_seg_len = len(get_valid_data(diff_data[month:cp]))
                            #left_seg_len = len(diff_data[month:cp])
                            print ("CHG1: ",neighb_index,neighb_id,"num,edit,1b,1e,imo,nhits",
                                   left_seg_len,"--",cp+1,left,month,neighb_hits)
                            
                            left = month
                            break
                    # ... and the final left-segment is too short
                    if left_seg_len < iminlen:
                        print ("Low1: ",neighb_index,neighb_id,"num,edit,1b,1e,imo,nhits",
                                left_seg_len,"--",cp+1,left,month,"--")
                        continue
                        
                    ## We can now estimate the raw changepoint amplitude using minbic.
                    ## However, we'll short-circuit a lot of the work by telling it to only 
                    ## use the KTHTPR0 model (simple step-change model)
                    (seg_x, seg_data) = range(left+1, right+1), diff_data[left+1:right+1]
                    bp_index = cp-(left+1)
                    #print left, cp, right, "|", bp_index
                    #print left_seg_len, right_seg_len
                    bic_result = minbic(seg_x, seg_data, bp_index, missing_val,
                                        models=[('KTHTPR0', kthtpr0), ])
                    ## Also check the first difference correlations between the 
                    ## monthly anomalies
                    station_first_diff = compute_first_diff(station_anomalies, missing_val)
                    neighb_first_diff = compute_first_diff(neighb_anomalies, missing_val)
                    corr = compute_corr(station_anomalies, neighb_anomalies)
                    
                    ## Write out the results of this testing process so far
                    cmodel = bic_result['cmodel']
                    bic = bic_result['bic']
                    test_stat = bic_result['test_stat']
                    crit_val = bic_result['crit_val']
                    offset = bic_result['offset']
                    slopes = bic_result['slopes']
                    left_slope, right_slope = slopes
                    print ("%s %6s-%6s %s %7.2f %7.2f %7.2f %7.2f %7.3f %7.3f -- %d --" % 
                           (tstr,id,neighb_id,chgptstr,crit_val,test_stat,offset,corr,
                            left_slope,right_slope,right_seg_len))
                    
                    ## Analysis is done.
                    ## Keep the adjustment (offset) for each neighbor/segment,
                    ## set/reset trend for each neighbor/segment
                    ##     the first segment is the left-segment,
                    ##     the second segment is the right-segment
                    ##
                    ## Note that we reset left/right potentially to avoid conflicts
                    ## within the paired neighbor data. However, our estimates of
                    ## trends/offsets associated with the "right" adjacent changepoint
                    ## actually refers to that original right changepoint. We'll 
                    ## reset left, cp, and right from the bracket before continuing
                    (left, cp, right) = bracket[:]
                    # Do the left segment first
                    left_dict = dict()
                    left_dict['adj'] = offset*correction
                    left_dict['cor'] = corr
                    left_dict['bic'] = bic
                    left_dict['cmodel'] = cmodel
                    left_dict['trend'] = left_slope
                    left_dict['spanob'] = left_seg_len
                    pair_results[neighb_id][cp] = left_dict
                                    
                    # Do the right segment now
                    right_dict = dict()
                    right_dict['adj'] = offset*correction
                    right_dict['cor'] = corr
                    right_dict['bic'] = bic
                    right_dict['cmodel'] = cmodel
                    right_dict['trend'] = right_slope
                    right_dict['spanob'] = right_seg_len
                    
                    if right not in pair_results[neighb_id]:
                        pair_results[neighb_id][right] = right_dict
                    else:
                        # We've already recorded this segment before for the last
                        # changepoint. Update the slopes/spanob count (length of
                        # preceding segment) if the slopes are different and the 
                        # length is different.
                        new_trend = slopes[1]
                        new_spanob = right_seg_len
                        old_trend = pair_results[neighb_id][right]['trend']
                        old_spanob = pair_results[neighb_id][right]['spanob']
                        if old_trend != new_trend:
                            print (" Seg2 diff: %s %4d old: %7.2f %4d new: %7.2f %4d" % 
                                   (pair_str,right,old_trend,old_spanob,new_trend,new_spanob))
                            # if the new count is greater than the old one, the slope
                            # is probably more robust so update those entries.
                            if new_spanob > old_spanob:
                                pair_results[neighb_id][right]['trend'] = new_trend
                                pair_results[neighb_id][right]['spanob'] = new_spanob
                                
                    ## We're done with this pair/changepoint. Summary output -
                    if step == 2:
                        print "itarg,ipair,ichg,numc,iqt,adj,trends: -- -- -- --",cmodel,offset,slopes
                #raw_input("pause")
                
                ####################################################################
                ## ADJUSTMENT DETERMINATION SECTION
                # Recall the paired-changepoint analyses we just performed, and 
                # determine if the potential adjustment is statistically valid
                (left, cp, right) = bracket[:]
                
                pair_data = []
                for neighb_id in pair_results:
                    if not cp in pair_results[neighb_id]:
                        continue
                    
                    cp_stats = pair_results[neighb_id][cp]
                    adjacent_stats = pair_results[neighb_id][right]
                    
                    trends = (cp_stats['trend'], adjacent_stats['trend'])
                    
                    pair_dict = dict(neighb_id=neighb_id, adj=cp_stats['adj'],
                                     cor=cp_stats['cor'], trends=trends, used=True)
                    pair_data.append(pair_dict)
                    
                npairs = len(pair_data)
                if npairs < numclim:
                    print "Adjpass numc low --",station_index,id,left,cp,right,npairs
                    continue
                
                # Process -
                #    1) Remove both adjustment and trend outliers
                #    2) Calculate median adjustment
                #    
                #    filter around inter-quartile range
                qscale = hom_params['qscale']
                pair_data = sorted(pair_data, key=operator.itemgetter('adj'))
                pair_chgs = [p['adj'] for p in pair_data]
                chg_25th, chg_median, chg_75th = tukey_med(pair_chgs)
                
                chg_iqr = chg_75th-chg_25th
                chg_low = chg_25th - (chg_median-chg_25th)*1.0*qscale
                chg_high = chg_75th + (chg_75th-chg_median)*1.0*qscale
                print (" TRIM p25, p75, pct50, rng, lo, hi: %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f" % 
                       (chg_25th,chg_75th,chg_median,chg_iqr,chg_low,chg_high) )         
                # If any of the estimated changepoints are outside the statistically
                # robust range we just computed, then flag them as we print them and
                for data in pair_data:
                    neighb_id = data['neighb_id']
                    neighb_index = all_station_list.index(neighb_id)
                    adj = data['adj']
                    cor = data['cor']
                    trends = data['trends']
                    
                    if not (chg_low < adj < chg_high):
                        data['used'] = False                
                    flag = 'U' if data['used'] else 'X'
                    print ("%s %4d %7.2f %8.4f %8.4f %7.2f" % 
                           (flag,neighb_index,adj,trends[0],trends[1],cor) )
                
                valid_adj_count = len([d for d in pair_data if d['used']])
                if valid_adj_count < numclim:
                    if step == 2:
                        print ("Insuff trimmed mean -- %4d %s %5d %5d %5d %5d" % 
                               (station_index,id,left,cp,right,valid_adj_count) )
                        continue
                
                ## BUG: The code here re-computes the inter-quartile range by 
                ##     scaling qscale by 1.0. Curiously, it doesn't reject any
                ##     pairs based on this new range.
                chg_iqr = chg_75th-chg_25th
                chg_low = chg_25th - (chg_median-chg_25th)*qscale
                chg_high = chg_75th + (chg_75th-chg_median)*qscale
                
                ## Tweak the inter-quartile range to check if the adjustment is 
                ## Check whether the computed adjustment is significant. That is,
                ## if 0 is included within the inter-quartile range we computed, then
                ## we can't reject the null hypothesis that the changepoint is significant
                if chg_high*chg_low > 0.0:
                    # signs are the same, so 0 isn't included in the range.
                    procstr = "CONSHF"
                    sigadj = chg_median
                else:
                    procstr = "ZERSHF"
                    sigadj = 0.0
                
                final_results[id][cp] = dict(adj=sigadj, std=chg_iqr*1.0*qscale,
                                             num=npairs)
    
                print ("%2d %s-%s %s %7.2f" % 
                       (station_index,id,procstr,chgptstr,sigadj) )
            
            ## Print some final output about what changepoints remain for this station
            final_station_results = final_results[id]
            final_cps = sorted(final_station_results.keys())
            for cp in final_cps:
                adj = final_station_results[cp]['adj']
                std = final_station_results[cp]['std']
                
                cp_stats = station_cp_dict[cp]
                hit_count = cp_stats['jsum']    
                iy, im = imo2iym(cp)
                
                print ("-- %5d %s Estamt chgout: -- %4d%2d %5d %5d %7.2f %7.2f" % 
                       (station_index+1,id,iy,im,cp,hit_count,adj,std) )
            #raw_input("pause")
            
        ## Remove the accumulated non-significant changepoints (either non-sig because
        ## there was too much missing data, the target segment was too short, or the 
        ## trimmed mean test could not reject the null hypothesis of no change
        for id in station_list:
            station_index = station_list.index(id)
            
            final_station_results = final_results[id]
            final_cps = sorted(final_station_results.keys())
            for cp in final_cps:
                iy, im = imo2iym(cp)
                cp_index = final_cps.index(cp)
                adj = final_station_results[cp]['adj']
                std = final_station_results[cp]['std']
                if adj == 0.0:
                    print ("%s %5d Remove chgpt %5d %4d %2d %4d" % 
                           (id,station_index,cp_index,iy,im,cp))
                    del network.raw_series[id].changepoints[cp]
                else: 
                    # Update the network's record of changepoints with this new list
                    network.raw_series[id].changepoints[cp]['ahigh'] = adj
                    network.raw_series[id].changepoints[cp]['astd'] = std
            # the changepoint at first month has been removed; add it back in
            network.raw_series[id].changepoints[first] = dict(ahigh=0.0,
                                                              astd=0.0,
                                                              jsum=0)   
    ## done!                       
                
def tukey_med(data):
    """Computes the 25th, 50th, and 75th percentiles of a set of data, using
    Tukey's median method. Assume the data is already sorted.
    
    """
    
    num_data = len(data)
    
    # Compute median (50th percentile)
    median_index = num_data/2
    if not num_data%2: # if odd...
        median = data[median_index]
    else:
        median = 0.5*(data[median_index-1]+data[median_index])
    
    # Compute lower-quarter (25th percentile), inclusive of the median
    if not median_index%2: 
        low_quarter = data[median_index/2]
    else:
        low_index = median_index-(median_index/2)
        low_quarter = 0.5*(data[low_index-1]+data[low_index])
        
    # Compute upper-quarter (75th percentile), inclusive of the median
    if not median_index%2:
        high_quarter = data[median_index+(median_index/2)]
    else:
        high_index = median_index+(median_index/2)
        high_quarter = 0.5*(data[high_index-1]+data[high_index])
        
    return low_quarter, median, high_quarter
    
    
    
                
            
        
        
        