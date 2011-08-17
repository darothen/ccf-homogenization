#!/usr/bin/env python
#
# Copyright (C) 2011 Daniel Rothenberg.
# See Google Code project page for license, 
# http://code.google.com/p/ccf-homogenization/
"""Filters to parse most likely undocumented changepoints from those detected
during splitmerge.

"""
__docformat__ = "restructuredtext"

import numpy as np

from util import imo2iym
from util import scale_series, compute_monthly_avg_std

def summary_table(station_ids, hits, all_data, print_missing=False, 
                  missing_val=-9999, **hom_params):
    """Print out a summary of a specified working array in matrix form. The 
    matrix will contain each entry in hits, and if print_missing is flagged
    as True, it will check all_data and replace an entry in hits with a mark
    that the value is missing in the original dataset.
    
    :Param station_ids:
        A list of the station_ids corresponding to the first axis in both
        hits and all_data.
    :Param hits:
        The working array to print in the matrix.
    :Param all_data:
        The array to check against for missing values.
    :Param print_missing:
        (optional) Enable printing of locations where "missing values" are
        detected in all_data, as denoted by the given missing_val
    :Param missing_val:
        (optional) Placeholder for missing values in all_data. Default is
        -9999.
    :Return:
        Nothing; prints the 'hits' array to console.
    
    """

    ## Print header - 
    head1 = "              |"+"|".join([i[:3] for i in station_ids])+"|"
    head2 = "              |"+"|".join([i[3:] for i in station_ids])+"|"
    print head1
    print head2
    
    ## Print monthly series basic
    def con2str(data, missing_val=-9999):
        """Convert entry in the working array into a string to print"""
        val, hits = data
        
        if val == missing_val and print_missing:
            return "-X-"
        elif hits > 0:
            return "%3d" % hits
        else:
            return "---"
        
    num_stations, num_months = all_data.shape
    for imo in xrange(num_months):
        year, month = imo2iym(imo)
        base_str = " %4d %2d %4d |" % (year, month, imo-11)
        month_strs = "|".join(map(con2str, 
                                  zip(all_data[:,imo],hits[:,imo])))+"|"
        print_month_strs = False
        for i in range(10):
            if str(i) in month_strs: 
                print_month_strs = True
                break
        
        if print_month_strs:
            print base_str+month_strs

def filter1(network, **hom_params):
    """Attempts to 'unconfound' the undocumented changepoints detected for each
    pair of stations and attribute them to a single station. This first filter
    also takes care of book-keeping by setting up arrays to hold the results of
    the filters employed in the Menne/Williams pairwise homogenization algorithm.
    
    :Param network:
        The network object containing the station/series data and metadata, as
        well as the paired results from the splitmerge process, stored in a 
        dictionary called pair_results,
        :Ivar pair_results:
            A dictionary containing the paired results. Should have keys of the
            form "%s-%s" % (id1, id2), and each key should be associated with a
            dictionary element. That dictionary has keys corresponding to the
            breakpoint and the associated analysis computed by minibic() for that
            breakpoint, as well as the keys "nspan" and "del", which are lists of
            months flagged as suspect during the splitmerge process.
    :Param hom_params:
        The parameters object containing important parameter settings for this
        analysis:
        :Ivar nmo:
            The number of months total in all of the raw series in the network.
    :Return:
        Updates the network object with the following fields:
        :Ivar hits_array:
            A numpy array of dimensions [number of stations, number of months]
            recording the filtered set of changepoints. The array will be 0 
            everywhere except where a station logs a changepoint, and in those
            instances, the array will reflect the number of times that 
            changepoint was implicated in a paired neighbor.
        :Ivar amps_array:
            A numpy array of dimensions [number of stations, number of months]
            recording the normalized changepoint amplitude estimated at a given
            month for a given station. If there is no breakpoint at a given
            entry, records 0.0; else, this is the average of the offsets
            for each paired hitpoint divided by their standard deviation.
        
    """

    ## Setup the data in NumPy arrays for easy inspection.
    station_list = network.stations.keys()
    ids = station_list
    all_data = [network.raw_series[id].monthly_series for id in ids]
    all_data = np.array(all_data)
        
    ################################################################################
    ## PRE FILTER 1
    ## Search through the record of pair results and initially attribute every
    ## detected changepoint to *both* stations in the pair. Also, delete spans
    ## of data which were flagged as suspect or troublesome during the changepoint
    ## detection process.    
    hits = np.zeros_like(all_data)
    
    hits_neighbors = [[list() for i in xrange(hom_params['nmo'])] for i in xrange(len(ids))]
    deleted = np.zeros_like(all_data)
    hits_models = np.zeros(all_data.shape, dtype=np.dtype( (str, 7) ))
    
    for pair in network.pair_results:
        id1, id2 = pair.split("-")
        id1_ind, id2_ind = ids.index(id1), ids.index(id2)
    
        result = network.pair_results[pair]
        for (bp, result) in result.iteritems():
            if bp == "nspan":
                continue
            elif bp == 'del':
                #continue
                ## Store in network
                network.raw_series[id1].delete_months(network.pair_results[pair]['del'])
                network.raw_series[id2].delete_months(network.pair_results[pair]['del'])            
                
                print "DEL", pair, network.pair_results[pair]['del']
                for bad_month in network.pair_results[pair]['del']:
                #    hits[id1_ind, bad_month] += 1
                #    hits[id2_ind, bad_month] += 1
                #    deleted[id1_ind, bad_month] += 1
                #    deleted[id2_ind, bad_month] += 1
                    pass
                
            else:
                hits[id1_ind, bp] += 1
                hits[id2_ind, bp] += 1
           
                hits_neighbors[id1_ind][bp].append(id2)
                hits_neighbors[id2_ind][bp].append(id1)
                
                hits_models[id1_ind, bp] = result['cmodel']
                hits_models[id2_ind, bp] = result['cmodel']

    ## Print a summary of the initial list of suspect changepoints
    summary_table(station_list, hits, all_data, False)    
    
    ################################################################################
    ## FILTER 1
    ## Unconfound the paired changepoints by repeatedly finding the most common 
    ## changepoint, and attributing that as the culprit for all of its pairs in 
    ## the month it occurs while deleting the hits the culprit has caused. Repeat
    ## until there is no changepoint which has been a culprit less than two times.
    new_hits = np.zeros_like(hits)
    amps = np.zeros_like(hits)
    
    ## Loop forward through all of the months, analyzing one at a time
    for imo in range(hits.shape[1]):
        
        hits_month = hits[:,imo]
        while hits_month.max() > 1:
            
            max_in_hits = hits_month.max()
            station_index = hits_month.argmax()
        
            station_id = station_list[station_index]
            iy, im = imo2iym(imo)

            ## Find entries in pair_results with this station and a changepoint on this
            ## date
            pr_keys = [key for key in network.pair_results if (station_id in key and
                                                               imo in network.pair_results[key])]

            if not pr_keys: 
                break
            
            offset_sum, offset_z_sum, count = 0.0, 0.0, 0.0
            for key in pr_keys:
                bp_summary = network.pair_results[key][imo]
    
                id1, id2 = key.split("-")
                
                offset, offset_z = bp_summary['offset'], bp_summary['offset_z']
                if station_id == id2: 
                    offset = offset*-1.0            
                    
                offset_sum += offset
                offset_z_sum += abs(offset_z)
                count += 1
            avg_offset = offset_sum/count
            avg_offset_z = offset_z_sum/count
            
            print ("  -- %s-CONFRM MW1 at %5d %4d %2d AVG ADJ: %3.2f %2.2f %3d" % 
                   (station_id, imo, iy, im, avg_offset, avg_offset_z, max_in_hits))
            
            new_hits[station_index,imo] = max_in_hits
            amps[station_index,imo] = avg_offset_z
            hits[station_index,imo] = 1
            hits_month[station_index] = 1
                    
    ## Print a summary of the filtered changepoints.
    summary_table(station_list, new_hits, all_data, False)
    
    ## Update network with this new data
    network.hits_array = new_hits
    network.amps_array = amps
    
def filter2(network, **hom_params):
    """Reconciles the detected undocumented changepoints with documented ones,
    if available, by "absorbing" detected changepoints near these known breaks
    in the data. 
    
    TODO: implement this functionality
    
    :Param network:
        The network object containing the station/series data and metadata, and
        a record of suspect, undocumented changepoints and their approximate
        amplitude:
        :Ivar hits_array:
            A numpy array of dimensions [number of stations, number of months]
            containing the timing of the set of undocumented changepoints 
            filtered so far. 
        :Ivar amps_array:
            A numpy array of dimensions [number of stations, number of months]
            containing the approximate amplitudes of the undocumented 
            changepoints filtered so far.
    :Param hom_params:
        The parameters object containing important parameter settings for this
        analysis:
        <none so far>
    :Return:
        does nothing right now other than print out some summary information.
        
    """
    
    ################################################################################
    ## FILTER 2
    ## won't actually do anything but print things to console for now
    ## Use station history and metadata to absorb undocumented changepoints
    ## to known ones where possible
    station_list = network.stations.keys()
    ids = station_list
    
    (num_stations, num_months) = network.amps_array.shape
    for station_index in range(num_stations):
        
        station_id = station_list[station_index]
        station_series = network.raw_series[station_id]
        data = station_series.series
        scale_series(data, 0.1, station_series.MISSING_VAL)
        stdk = compute_monthly_avg_std(data)
        
        for month in range(num_months):
            
            if network.hits_array[station_index, month] > 0:
                ahigh = network.amps_array[station_index, month]
                jsum = network.hits_array[station_index, month]
                
                astd = abs(ahigh)
                
                iy, im = imo2iym(month)
                
                ## These are PRE-DEFINED in inhomog.parm.system.mthly.incl. They 
                ## will need to be re-factored later in to a more logical place
                ## (Parameters maybe?)
                arange = [0.4, 0.6, 0.8, 1.0, 1.5, 3.0, 5.0]
                mrgyr = [36, 18, 12, 8, 6, 5, 5]
                nrange = len(arange)
                
                # Create search bracket for looking at station history files
                for irange in range(len(arange)):
                    if astd < arange[irange]: break
                ibracket = mrgyr[irange]*2 + 1
                
                print "ASTD: ",station_index,station_id,"1",iy,im,astd,ahigh,stdk,jsum,ibracket
    
def filter3(network, **hom_params):
    """Reduce the number of detected changepoints by condensing closely-related
    changepoints into one another, by choosing the biggest amplitude change.
    
    :Param network:
        The network object containing the station/series data and metadata, and
        a record of suspect, undocumented changepoints and their approximate
        amplitude:
        :Ivar hits_array:
            A numpy array of dimensions [number of stations, number of months]
            containing the timing of the set of undocumented changepoints 
            filtered so far. 
        :Ivar amps_array:
            A numpy array of dimensions [number of stations, number of months]
            containing the approximate amplitudes of the undocumented 
            changepoints filtered so far.
    :Param hom_params:
        The parameters object containing important parameter settings for this
        analysis:
        <none so far>
    :Return:
        Record the final set of changepoints for each series in the given
        network object as the following field:
        :Ivar changepoints:
            A dictionary containing the changepoint locations and amplitudes. 
            Each key in the dictionary is the integer month index where a 
            breakpoint occurs, and each dictionary has the following fields:
            :Param jsum:
                The number of times this changepoint caused a break to be logged
                in a paired neighbor during splitmerge
            :Param ahigh:
                The amplitude change associated with the changepoint
            :Param astd:
                The error associated with the amplitude change
    """
    
    ################################################################################
    ## FILTER 3 TEST
    ## Try to reduce the number of changepoints by condensing closely-related
    ## changepoints into one another.
    ##
    ## Go back through the years/months
    ##    For each technique
    ##        find the highest remaining hits
    ##        accumulate all of the new_hits and ntests for each month +/- mrgyr 
    ##           (from amps) while skipping missing data
    ##    All new_hits and amps are zeroed out
    ##    For all accumulations within nmrgryr
    ##        when given month are greater or equal to ithresh then add to
    ##           the nhits and ntests arrays
    ##
    ## iconfirm := 2 | the min number of hits for a possible breakpoint to be valid
    ## nmrgyr := -2 | no idea what this does; it defaults to -2 
    station_list = network.stations.keys()
    ids = station_list
        
    inconfirm = 2
    nmrgyr = -2
    
    final_hits = np.zeros_like(network.hits_array)
    ifound = 0

    (num_stations, num_months) = network.amps_array.shape
    for station_index in range(num_stations):
        
        station_id = station_list[station_index]
        data = network.raw_series[station_id].series
        miss = network.raw_series[station_id].MISSING_VAL
        data_monthly = network.raw_series[station_id].monthly_series
        stdk = compute_monthly_avg_std(data)
        
        ## setup temp arrays 
        khits = np.zeros(hom_params['nmo'])
        ktests = np.zeros(hom_params['nmo'])
        akhigh = np.zeros(hom_params['nmo'])
        
        # iterate until there are no more high points
        istop = False
        while not istop:
            
            ihighit, ahigh = 0.0, 0.0
            
            # find the highest count - the most number of hits at a possible breakpoint
            for month in range(num_months):
                
                isum, asum = 0.0, 0.0
                if network.hits_array[station_index, month] >= inconfirm:
                    jhit = network.hits_array[station_index, month]
                    isum += jhit
                    asum += network.amps_array[station_index, month]*jhit
                
                # find the highest chgpt hit station by hits
                if isum > ihighit: 
                    ihighit = isum
                    ihighmo = month
                    ahigh = asum / isum 
            ## now -
            ##    ihighmo := month of highest hit value
            ##    ihighit := sum of hits over all tests
            ##    ahigh := estimated adjustment
            print "----itarg,ihighit,ihighmo,ahigh,stdk",station_index,ihighit,ihighmo,ahigh,stdk
                    
            #Keep going until there are no more hits
            if ihighit > 0:
                # bracket the highest hit +/- nmrgyr
                if nmrgyr != -2:
                    ibracket = nmrgyr*2 + 1
                else:
                    # else bracket using amplitude of chgpt to define month range
                    ## This are PRE-DEFINED in inhomog.parm.system.mthly.incl. They 
                    ## will need to be re-factored later in to a more logical place
                    ## (Parameters maybe?)
                    arange = [0.4, 0.6, 0.8, 1.0, 1.5, 3.0, 5.0]
                    mrgyr = [36, 18, 12, 8, 6, 5, 5]
                    nrange = len(arange)
                    
                    # Create search bracket for looking at station history files
                    for irange in range(len(arange)):
                        astd = abs(ihighit)
                        if astd < arange[irange]: break
                    ibracket = mrgyr[irange]*2 + 1       
                
                # go through the bracket, look for the highest hits already found
                # start at the hit-point index and expand outward in the series
                # i.e, if our ihighmo = 10, look at [10, 11, 9, 12, 8, 13, 7, ...]
                # keep track of missing values in both directions.
                max_radius = ibracket/2
                miss_left, miss_right = 0, 0
                radius = 0
                absorbed = False
                while radius < max_radius:
                    
                    ## DEAL WITH THE RIGHT MONTH
                    right_month = ihighmo + radius + miss_right
                    if right_month == hom_params['nmo']:
                        break
                    if data_monthly[right_month] == miss:
                        # this month is missing, go to the next
                        while data_monthly[right_month] == miss:
                            right_month += 1
                    # Absorb lesser hit into the closest higher hit
                    #print right_month, right_month-5, right_month
                    #print np.where(khits > 0)
                    if khits[right_month] > 0:
                        khits[right_month] += ihighit
                        akhigh[right_month] += ahigh*ihighit
                        print "Absorb hit: ",station_index,ihighmo," to ",right_month,khits[right_month],ktests[right_month],akhigh[right_month]/khits[right_month]
                     
                        # zero test array block for next iter
                        network.hits_array[station_index,ihighmo] = 0
                        network.amps_array[station_index,ihighmo] = 0.0
                        absorbed = True
                        break
                     
                    ## DEAL WITH THE LEFT MONTH
                    left_month = ihighmo - radius - miss_left
                    if left_month == 0:
                        break
                    if data_monthly[left_month] == miss:
                        # this month is missing, go to the next
                        while data_monthly[left_month] == miss:
                            left_month -= 1
                    
                    # Absorb lesser hit into the closest higher hit
                    #print left_month, left_month-5, left_month+5
                    #print np.where(khits > 0)
                    if khits[left_month] > 0:
                        khits[left_month] += ihighit
                        akhigh[left_month] += ahigh*ihighit
                        print "Absorb hit: ",station_index,ihighmo," to ",left_month,khits[left_month],ktests[right_month],akhigh[right_month]/khits[right_month]
                        
                        # zero test array block for next iter
                        network.hits_array[station_index,ihighmo] = 0
                        network.amps_array[station_index,ihighmo] = 0.0
                        absorbed = True
                        break
                
                    radius += 1
                        
                # if no hits found, setup new hit
                if not absorbed:
                    khits[ihighmo] = ihighit
                    ktests[ihighmo] = 1
                    akhigh[ihighmo] = ahigh*ihighit
                    print "New CHG hit: ",station_index,ihighmo,khits[ihighmo],ktests[ihighmo],akhigh[ihighmo]/khits[ihighmo]
            
                    network.hits_array[station_index, ihighmo] = 0
                    network.amps_array[station_index, ihighmo] = 0.0
                
                #raw_input("pause")
            else:
                istop = True
                
        print "----------------------------------------------"
        # examine interim khits array for station's filtered changepoints
        uchgpt_dict = dict()
        for month in range(hom_params['nmo']):
            # ... if highest hits > ithres(npair) then save
            # fetch the numbr of pairs tested
            
            if khits[month] > 0:
                npair = ktests[month]
                jsum = khits[month]
                
                ihthres = 2
                iy,im = imo2iym(month)
                print "itarg,imo,iym,npair,jsum,ihthres,stdk",station_index,month,iy,im,npair,jsum,ihthres,stdk
                if jsum >= ihthres:
                    # passed threshold test- put interim into final
                    final_hits[station_index, month] += jsum
                    ifound += 1
                    
                    # debug stuff
                    ahigh = akhigh[month]/khits[month]
                    astd = ahigh
                    print ("%5d %6s-UCHGPT KW%1d at %4d %3d %4d %6.2f %6.2f %3d %3d %3d" % 
                           (station_index, station_id, 1, iy, im, jsum, ahigh, astd, ibracket, npair, ihthres) )
                    
                    uchgpt_dict[month] = { 'jsum': jsum,
                                           'ahigh': ahigh,
                                           'astd': astd }
                    
        network.raw_series[station_id].changepoints = uchgpt_dict
        
    print "-------------------------------------------------"
    print "Undoc filter: ",ifound
    
""" Test code for Filter 1 - ignore
##istop = False
##while not istop:
##    ## Find the highest chgpt occurrence for all targets/all yr-mths
##    high_hits = 0
##    for month in range(hom_params.nmo):
##        for id_ind in range(hom_params.nstns):
##            if hits[id_ind, month] > high_hits:
##                high_hits = hits[id_ind, month]
##                high_station = id_ind
##                high_month = month
##                
###    high_hits = hits.max()
###    flat_index = hits.argmax()
###    high_station, high_month = np.unravel_index(flat_index, hits.shape)
##    high_station_id = station_list[high_station]
##    
##    if high_hits > 1:
##        iy, im = imo2iym(high_month)
##        print " Yr/mth: ",high_month,iy,im," ihighnet: ",high_station," ihighit: ",high_hits
##        
##        for pair_id in station_list:
##            ida, idb = sorted([high_station_id, pair_id])
##            pair_key = "%s-%s" % (ida, idb)
##            pair_index = station_list.index(pair_id)
##            
##            if pair_key not in pair_results:
##                print " No link from ",high_station_id," back to ",pair_id
##                continue
##                    
##            if high_month in pair_results[pair_key]['del']:
##                continue          
##
##            pair_nspan = pair_results[pair_key]['nspan']
##            nchg = pair_nspan[high_month]
##            # Go backward... clearing out peripherals
##            for month in range(high_month-1,0,-1):
##                print "bkt imo: ", month
##                if nchg != pair_nspan[month]:
##                    break
##                print "      nspan: ",pair_results[pair_key]['nspan'][month]," ifndshow: ",hits[high_station,month]
##                pair_results[pair_key]['nspan'][month] = 0
##                hits[high_station, month] -= 1
##            month1 = month
##            # Go forward... clearing out peripherals
##            for month in range(high_month, hom_params.nmo):
##                print "fwt imo: ", month
##                if nchg != pair_nspan[month]:
##                    break
##                print "      nspan: ",pair_results[pair_key]['nspan'][month]," ifndshow: ",hits[high_station,month]
##                pair_results[pair_key]['nspan'][month] = 0
##                hits[high_station, month] -= 1  
##            month2 = month
##        #######################################################################
##        ## Repeat the clearing peripherals process for the neighbor station
##            # Go backward... clearing out peripherals
##            for month in range(high_month,0,-1):
##                print "bkp imo: ",month
##                if nchg != pair_nspan[month]:
##                    #print month, pair_results[pair_key].keys()
##                    break
##                print "      nspan: ",pair_results[pair_key]['nspan'][month]," ifndshow: ",hits[high_station,month]
##                pair_results[pair_key]['nspan'][month] = 0
##                hits[pair_index, month] -= 1
##            month1 = month
##            # Go forward (including the high month)... clearing out peripherals
##            for month in range(high_month, hom_params.nmo):
##                print "fwp imo: ",month
##                if nchg != pair_nspan[month]:
##                    #print month, pair_results[pair_key].keys()
##                    break
##                print "      nspan: ",pair_results[pair_key]['nspan'][month]," ifndshow: ",hits[high_station,month]
##                pair_results[pair_key]['nspan'][month] = 0
##                hits[pair_index, month] -= 1  
##            month2 = month
##            print "Zeroed STN: ",high_station_id," neigh:",pair_id," neigh index:",pair_index," range:",month1,month2," fndshow:",hits[pair_index,month]
##        new_hits[high_station, high_month] = hits[high_station, high_month]
##        hits[high_station, high_month] *= -1
##    else:
##        istop = True
#### All of the peripheral hits associated with each filtered changepoint should
#### be zeroed out now, leaving only attributed chgpt-station data.
#### Go through years/months again...
##for month in range(hom_params.nmo):
##    for station_index in range(hom_params.nstns):
##        
##        if new_hits[station_index, month] > 1:
##            local_max = new_hits[station_index, month]
##            station_id = station_list[station_index]
##            
##            pr_keys = [key for key in pair_results if (station_id in key and
##                                                       month in pair_results[key])]
##            
##            if not pr_keys: break
##            
##            offset_sum, offset_z_sum, count = 0.0, 0.0, 0.0
##            for key in pr_keys:
##                bp_summary = pair_results[key][month]
##    
##                id1, id2 = key.split("-")
##                
##                offset, offset_z = bp_summary['offset'], bp_summary['offset_z']
##                if station_id == id2: 
##                    offset = offset*-1.0            
##                    
##                offset_sum += offset
##                offset_z_sum += abs(offset_z)
##                count += 1
##            avg_offset = offset_sum/count
##            avg_offset_z = offset_z_sum/count
##            amps[station_index, month] = avg_offset_z
##            
##            iy, im = imo2iym(month)
##            print ("  -- %s-CONFRM MW1 at %5d %4d %2d AVG ADJ: %3.2f %2.2f %3d" % 
##                   (station_id, month, iy, im, avg_offset, avg_offset_z, local_max))
"""
    
    