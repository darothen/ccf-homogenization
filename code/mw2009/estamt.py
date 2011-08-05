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
    
    station_list = network.stations.keys()
    all_data = [network.raw_series[id].monthly_series for id in station_list]
    all_data = np.array(all_data)
    
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
    ## STEP 2 - NOT SIG REMOVAL
    ## equivalent to ipass loopback for istep == 2 in Fortran PHA
    print " ---------------- NOT SIG REMOVAL --------------- "
    tstr = "Not sig: "
    outid = "NS"
    ipass = 1
    iminlen = hom_params['minlen']
    numclim = 3
    
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
        for other_id in station_list:
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
            
        ## STATION CHGPT LOOP
        # Loop over "brackets" of changepoints - that is, for changepoints
        # [a, b, c, d], consider the two brackets [a,b,c] and [b,c,d] with 
        # the center value of the changepoints. Note that in the Fortran PHA,
        # we go through these brackets in reverse order - right to left.
        brackets = zip(sorted_cps[-3::-1], sorted_cps[-2::-1], sorted_cps[::-1])
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
                if id == id1:
                    neighb_id = id2
                else:
                    neighb_id = id1
                
                # Add this pair to pair_results if it's not already there
                (ida, idb) = sorted([id1,id2])
                pair_str = "%s-%s" % (ida, idb)
                if pair_str not in pair_results:
                    pair_results[pair_str] = dict()
                print pair_str
                
                neighb_index = station_list.index(neighb_id)
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
                # Do the left segment first
                left_dict = dict()
                left_dict['adj'] = offset
                left_dict['cor'] = corr
                left_dict['bic'] = bic
                left_dict['cmodel'] = cmodel
                left_dict['trend'] = left_slope
                left_dict['spanob'] = left_seg_len
                pair_results[pair_str][cp] = left_dict
                                
                # Do the right segment now
                right_dict = dict()
                right_dict['adj'] = offset
                right_dict['cor'] = corr
                right_dict['bic'] = bic
                right_dict['cmodel'] = cmodel
                right_dict['trend'] = right_slope
                right_dict['spanob'] = right_seg_len
                
                if right not in pair_results[pair_str]:
                    pair_results[pair_str][right] = right_dict
                else:
                    # We've already recorded this segment before for the last
                    # changepoint. Update the slopes/spanob count (length of
                    # preceding segment) if the slopes are different and the 
                    # length is different.
                    new_trend = slopes[1]
                    new_spanob = right_seg_len
                    old_trend = pair_results[pair_str][right]['trend']
                    old_spanob = pair_results[pair_str][right]['spanob']
                    if old_trend != new_trend:
                        print (" Seg2 diff: %s %4d old: %7.2f %4d new: %7.2f %4d" % 
                               (pair_str,right,old_trend,old_spanob,new_trend,new_spanob))
                        # if the new count is greater than the old one, the slope
                        # is probably more robust so update those entries.
                        if new_spanob > old_spanob:
                            pair_results[pair_str][right]['trend'] = new_trend
                            pair_results[pair_str][right]['spanob'] = new_spanob
                            
                ## We're done with this pair/changepoint. Summary output -
                print "itarg,ipair,ichg,numc,iqt,adj,trends: -- -- -- --",cmodel,offset,slopes
            #raw_input("pause")
            
                    
                
                
            
        
        
        