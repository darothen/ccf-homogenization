#!/usr/bin/env python
#
# Copyright (C) 2011 Daniel Rothenberg.
# See Google Code project page for license, 
# http://code.google.com/p/ccf-homogenization/

"""Final stage of MW2009 process - Estimation of adjusted series

"""
__docformat__ = "restructuredtext"

from util import imo2iym

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
                        del network.raw_series[id].changepoints[del_key]
                    
                    del_str = "Del 1st segment: " if cp1 == first else "Delete segment: "
                    
                    print id,station_index+1,del_str,imo2iym(cp1),cp1,imo2iym(cp2),cp2            
            
            new_changepoints = network.raw_series[id].changepoints
            new_cps = sorted(new_changepoints.keys())
            print "  First data value: ",imo2iym(first)
            for cp in new_cps:
                print "  End seg:",new_cps.index(cp)," ym: ",imo2iym(cp),cp,new_changepoints[cp]['jsum']
            print "    End segment ym: ",imo2iym(last),last
            ####################################################################
                
