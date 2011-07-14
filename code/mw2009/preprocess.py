#!/usr/bin/env python
#
# Daniel Rothenberg, 2011-06-28

from util import compute_arc_dist, compute_first_diff, compute_corr, compute_std

from operator import itemgetter

def find_neighborhood(station, stations_list, numsrt=40, mindist=200.0, **kwargs):
    """
    Determines the neighborhood comprised of the closes numsrt neighbors located
    around a given station. 
    
    :Param stations:
        The candidate *Station* object.
    :Param stations_list:
        A list of *Station* objects to use to determine the neighborhood.
    :Param numsrt:
        The number of *Station* objects which should comprise the resulting 
        neighberhood.
    :Param mindist: 
        The minimum distance to use while searching through the *Station* objects
        iteratively. Most likely will not be used.
    
    :Return:
        The neighborhood as a list of two-element tuples. The tuples will each
        be of the form (Station.coop_id, computed distance from candidate), and
        the list will be sorted in descending order based on the computed
        distance.
    """
    
    neighbors = dict()
    for other in stations_list:
        
        if not station.coop_id == other.coop_id:
            dist = compute_arc_dist(station, other)
            neighbors[other.coop_id] = dist
        
    # Sort the dictionary of neighbors into a list, sorting from least
    # distance to greatest distances. Only take numsrt neighbors maximum.
    sorted_neighbors = sorted(neighbors.iteritems(),
                              key=itemgetter(1))[:numsrt-1]
    '''
    # This code is actually superfluous, but it loops through the list of
    # sorted neighbors, and finds only the neighbors which are within a maximum
    # distance from the candidate station. When its done, it sees if there are
    # enough neighbors, and if not, ups the threshold for distance and tries 
    # again.
    close_neighbors = []
    iter = 1
    hidist = mindist
    print "...searching for neighbors within %4.1fkm" % mindist
    while (len(close_neighbors) < numsrt-1):
        print "......iteration %d, %d neighbors found" % (iter, 
                                                          len(close_neighbors))
        for (coop_id, dist) in sorted_neighbors:
            if dist < hidist:
                sorted_neighbors.remove((coop_id, dist))
                close_neighbors.append((coop_id, dist))
        hidist = hidist+distinc
        print ".........distance now %4.1fkm" % mindist    
        iter = iter+1
    print "......ultimately found %d neighbors" % len(close_neighbors)
    all_neighbors[coop_id1] = close_neighbors
    '''

    return sorted_neighbors

def neighborhood_strings(station, neighborhood, stations_list):
    """
    Creates the output strings necessary to create a distance neighborhood file
    similar to what is produced by ushcn_dist_2004.v3. For each station in this
    type of file, there are 3 formatted lines, containing (in order) the station
    ids of the this station and its neighbors, the Fortran pointer indices
    corresponding to the neighbors, and the distances from this station to each
    of those neighbors (with the distance from this station to itself  
    defined as 0.0).
    
    This will be useful as a diagnostic tool in checking the new PHA code.
    
    :Param station:
        The candidate *Station* object whose neighborhood is being printed.
    :Param neighborhood:
        The sorted list of neighbors for the candidate *Station*. This should
        be computed by find_neighborhood.
    :Param stations_list:
        The list of *Station* objects from which the neighborhood was pulled. We
        only need this for the pointers string.
    :Return:
        A 3-tuple with the output strings for this neighborhood - a string with
        the station ids, a string with the pointer indices, and a string with 
        the station distances.
    
    """
    
    ids, dists = zip(*neighborhood)
    
    id_str = "%6s " % station.coop_id
    id_str = id_str + "".join(("%6s " % id for id in ids)) + "\n"
    
    ptr_str = "{0: >6d} ".format(stations_list.index(station)+1)
    other_ptrs = range(len(ids))
    ptr_str = ptr_str + "".join(("{0: >6d} ".format(ptr) for ptr in other_ptrs)) + "\n"
    
    dist_str = "   0.0 "
    dist_str = dist_str + "".join(("{0: >6.1f} ".format(dist) for dist in dists)) + "\n"
    
    return (id_str, ptr_str, dist_str)

def find_correlations(cand_series, series_dict, neighborhood, 
                      corrlim=0.1, begyr=1900, endyr=2010, minpair=14, 
                      numcorr=20):
    
    corr_dict = dict()
    
    coop_id1 = cand_series.coop_id
    print "...%s" % coop_id1
    neighbors = [id for id, dist in neighborhood]
    
    for coop_id2 in neighbors:
        print "......%s" % coop_id2
        neighb_series = series_dict[coop_id2]
        
        # Get the data for these series. Note that up until this point,
        # we haven't corrected for the fact that temperature data is reported
        # in tenths of a degree in the USHCN database. Let's go ahead and
        # correct that factor; it turns out that if you don't, the correlation
        # doesn't work correctly. Note that computing anomalies is a linear
        # operation, so it doesn't matter for the math so far that we've used
        # tenths of a degree instead of whole degrees.
        cand_data = [val*.1 for val in cand_series.monthly_series]
        neighb_data = [val*.1 for val in neighb_series.monthly_series]
        
        # We SHOULD have read the same years of data, and have equal lengths
        # of data series.
        assert cand_series.years == neighb_series.years
        assert len(cand_series.series) == len(neighb_series.series)
        
        # What is the missing value placeholder? Correct for being in tenths
        # of a degree.
        MISS = cand_series.MISSING_VAL*.1
        
        # Align the candidate and network series by looping through every
        # value, and choosing only months where BOTH a candidate and neighbor
        # value are present. If either or both are missing, skip that month
        # and go on to the next.
        print ".........Aligning cand/neighb series"
        cand_align, neighb_align = [], []
        for (cand_val, neighb_val) in zip(cand_data, neighb_data):
            if (cand_val != MISS and neighb_val != MISS):
                cand_align.append(cand_val)
                neighb_align.append(neighb_val)
        assert len(cand_align) == len(neighb_align)
                
        # We perform the correlation test on a first-difference series, so
        # compute that now. See util.compute_first_diff() for information on
        # what this operation entails.
        print ".........Computing first differences"
        cand_dif = compute_first_diff(cand_align, MISS)
        neighb_dif = compute_first_diff(neighb_align, MISS)
        
        # Now, we can actually compute the correlation coefficient. Again,
        # see util.compute_corr() for info on this mathematical operation.
        print ".........Computing correlation coefficient"
        r = compute_corr(cand_dif, neighb_dif, MISS, aligned=True)
        #r = compute_corr(cand_align, neighb_align, MISS)
        cand_std = compute_std(cand_dif, MISS)
        neighb_std = compute_std(neighb_dif, MISS)
        
        # If the correlation is above a threshold, we will keep it. In the
        # ushcn_corr_2004.v3 code, this threshold is 0.10.
        if r:
            print "            %1.3f %3.3f %3.3f" % (r, cand_std, neighb_std)
            corr_dict[coop_id2] = r

        else:
            print "            poor or no correlation"
           
    sort_corrs = sorted(corr_dict.iteritems(),  
                        key=itemgetter(1), reverse=True)
    good_corrs = [coop_id2 for (coop_id2, r) in sort_corrs if r > corrlim]
           
    nmonths = (endyr-begyr)*12
    ksum = [0]*nmonths
    jsum = [0]*nmonths
    lowtoo = [0]*nmonths
    kstns = 0
    # Determine ksum[nmonths], the number of neighbor data available to use 
    # in homogenizing data for this station at each month
    for imo in xrange(nmonths):
        if cand_data[imo] != MISS:            
            for (k, coop_id2) in zip(xrange(len(good_corrs)), good_corrs):
                neighb_series = series_dict[coop_id2]
                neighb_data = neighb_series.monthly_series
                if neighb_data[imo] != neighb_series.MISSING_VAL: 
                    ksum[imo] = ksum[imo]+1
            kstns = k
            jsum[imo] = ksum[imo]*1
                        
            if ksum[imo] < minpair:
                print " Total less than minpair: ",coop_id1,1900+(imo/12),1+(imo%12)
                lowtoo[imo] = 1
                
    # If we have more neighbors than necessary, then let's see if we can adjust
    # the numbers somewhat to bolster the amount of data in low-info periods,
    # being careful not too delete other good data.
    useful_neighbors = [n for n in good_corrs]
    jstns = kstns*1
    if kstns > numcorr-1:
        
        good_corrs.reverse()
        for (k, coop_id2) in zip(xrange(len(good_corrs)), good_corrs):
            
            iremove = 1
            npair = 0
            neighb_series= series_dict[coop_id2]
            neighb_data = neighb_series.monthly_series
            
            imonths = xrange(nmonths)
            CMISS, NMISS = MISS, neighb_series.MISSING_VAL
            
            iter_head = zip(imonths, cand_data, neighb_data)
            for (imo, c, n) in [(imo, c, n) for (imo, c, n) in iter_head]:
                if (c != CMISS) and (n != NMISS):
                    npair = npair+1
                    if ksum[imo] <= minpair:
                        print " Cannot remove:", coop_id1,'-',coop_id2,1900+(imo/12),1+(imo%12),ksum[imo],lowtoo[imo]
                        iremove = 0
                        break
                
            if iremove == 1:
                if kstns >= numcorr-1:
                    print " Remove:",coop_id1,"-",coop_id2,npair,corr_dict[coop_id2]
                    kstns = kstns-1
                    useful_neighbors.remove(coop_id2)
                    for imo in xrange(nmonths):
                        if cand_data[imo]!=CMISS and neighb_data[imo]!=NMISS:
                            ksum[imo] = ksum[imo]-1
        
    for imo in xrange(nmonths):
        if jsum[imo]>0:
            print "Original-Final:",coop_id1,1900+(imo/12),1+(imo%12),jsum[imo],ksum[imo]
    
    print "Original-Final Number stns:",coop_id1,jstns,kstns
    
    # Now, we know which neighbors a) are highly correlated with this station,
    # and b) add information where it is scarce in the temperature record.
    for coop_id2 in corr_dict.keys():
        if not coop_id2 in useful_neighbors:
            del corr_dict[coop_id2]
    
    return corr_dict
    
def tuple_list_to_dict(tuple_list):
    """
    Converts a list of two-element tuples of the form (key, value) into a 
    dictionary containing the same data.
    """
    new_dict = dict()
    for tuple in tuple_list:
        (key, value) = tuple
        new_dict[key] = value
        
    return new_dict

