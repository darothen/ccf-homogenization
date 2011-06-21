#!/usr/bin/env python
#
# Daniel Rothenberg, 2011-06-08

"""Sample driver script to compute network distances and correlations 
for feeding into Fortran PHA code"""

# http://docs.python.org/library/pprint.html
import pprint
# http://docs.python.org/library/random.html
import random
# http://docs.python.org/library/operator.html
from operator import itemgetter

# ccf-homogenization imports
import ushcn_io
from util import compute_arc_dist, compute_monthly_anomalies
from util import compute_first_diff, compute_corr, compute_std
import parameters


# A sample network of 52 stations. Started with a selection of literally
# 52 random stations from all the USHCN Coop stations, and then added
# Cheesman and Chula Vista because the MW2009 paper actually plots those
# two stations' monthly anomalies so I can make eyeball comparison
test_stations = ['215887', '041912', '034572', '116738', '361354',
                 '331152', '314684', '244522', '030936', '042941',
                 '331890', '291515', '098535', '258133', '266779', 
                 '324418', '427260', '086997', '449151', '181750', 
                 '307633', '200230', '124837', '200146', '220488', 
                 '413183', '215615', '415618', '461330', '206300', 
                 '200779', '121747', '046399', '231037', '431243', 
                 '411000', '302129', '240364', '111280', '043875', 
                 '298107', '234825', '228374', '164407', '248597', 
                 '315177', '443192', '314055', '153430', '120177', 
                 '041758', '051528']

# A Parameters object which contains the configuration for the homogenization
# analysis - specifically, the station list, settings for determining close
# neighbors, and settings for finding optimally correlated neighbors.
default_params = dict(nstns=50, 
                      mindist=200.0, 
                      distinc=200.0, 
                      numsrt=40,
                      numcorr=20, 
                      begyr=1900,
                      endyr=2010, 
                      data_src="raw", 
                      variable="max", 
                      stations=test_stations,
                      corrlim=0.1,
                      minpair=14)
params = parameters.default_parameters(**default_params)

pprint.pprint(params)

# Read in station data (download if necessary)
all_series, all_stations = ushcn_io.get_ushcn_data(params)

if not params.stations:
    station_ids = sorted(random.sample(all_stations.keys(), params.nstns))
else:
    station_ids = sorted(params.stations)
stations = dict(zip(station_ids, [all_stations[s] for s in station_ids]))

series_list = [all_series[station] for station in station_ids]
series = dict(zip([s.coop_id for s in series_list], series_list))

##########################################################################
print "Analyzing geographic network neighborhoods"

all_neighbors = dict()
for coop_id1 in station_ids:
    
    print coop_id1
    print "...computing neighbor distances"
    
    neighbor_dict = dict()
    for coop_id2 in station_ids:
        
        if not coop_id1 == coop_id2:
            dist = compute_arc_dist(stations[coop_id1], stations[coop_id2])
            #print "%s-%s = %4.3fkm" % (coop_id1, coop_id2, dist) 
            neighbor_dict[coop_id2] = dist

    # Sort the dictionary of neighbors into a list, sorting from least
    # distance to greatest distances. Only take numsrt neighbors maximum.
    sorted_neighbors = sorted(neighbor_dict.iteritems(),
                              key=itemgetter(1))[:params.numsrt-1]
    
    # This code is actually superfluous, but it loops through the list of
    # sorted neighbors, and finds only the neighbors which are within a maximum
    # distance from the candidate station. When its done, it sees if there are
    # enough neighbors, and if not, ups the threshold for distance and tries 
    # again.
    close_neighbors = []
    iter = 1
    hidist = params.mindist
    print "...searching for neighbors within %4.1fkm" % params.mindist
    while (len(close_neighbors) < params.numsrt-1):
        print "......iteration %d, %d neighbors found" % (iter, 
                                                          len(close_neighbors))
        for (coop_id, dist) in sorted_neighbors:
            if dist < hidist:
                sorted_neighbors.remove((coop_id, dist))
                close_neighbors.append((coop_id, dist))
        hidist = hidist+params.distinc
        print ".........distance now %4.1fkm" % params.mindist    
        iter = iter+1
    print "......ultimately found %d neighbors" % len(close_neighbors)
    all_neighbors[coop_id1] = close_neighbors
    
# Write a neighborhood output file. Since my algorithm for computing distance
# is slightly different than ushcn_dist_2004.v3, it won't produce exactly the
# same distance output file. However, all the distances are within 10km, which
# is perfectly fine. A bigger problem is that the ushcn_dist_2004.v3 outputs a
# list of pointers for referencing the various stations. Since my code doesn't
# emulate the Fortran code in terms of having arrays with lengths hard-coded,
# I don't pass around pointers in the same way, so the ptr_str produced and
# written here is garbage. It should not be a problem for coding MW2009, though,
# because I can find information about stations dynamically and easily.
print "...Assembling neighborhood output file"
dist_out = open("dist_out", 'wb')
for sta_id in test_stations:
    sorted_neighbors = sorted(all_neighbors[sta_id],
                              key=itemgetter(1))
    ids, dists = zip(*sorted_neighbors)
    id_str = ("%6s " % sta_id)+"".join(("%6s " % id for id in ids))+"\n"
    ptr_str = ("{0: >6d} ".format(test_stations.index(sta_id)+1))+"".join(("{0: >6d} ".format(test_stations.index(id)+1) for id in ids))+"\n"
    dist_str = ("   0.0 ")+"".join(("{0: >6.1f} ".format(d) for d in dists))+"\n"
    
    dist_out.writelines((id_str, ptr_str, dist_str))
dist_out.close()
    
##########################################################################

# Go through all the data we have, and replace the read-in values with
# monthly anomalies. Then, flatten the data into a list with all the data
# and length (endyr-begyr)*12
for s in series.itervalues():
    data = s.series
    anomalies = compute_monthly_anomalies(data, -9999)
    s.set_series(anomalies, s.years)
    
print "Determining correlated neighbors"

all_corrs = dict()
for coop_id1 in station_ids:
    
    corr_dict = dict()
    
    print "...%s" % coop_id1
    neighbors = [n[0] for n in all_neighbors[coop_id1]]
    for coop_id2 in neighbors:
        print "......%s" % coop_id2
        
        # Get the data for these series. Note that up until this point,
        # we haven't corrected for the fact that temperature data is reported
        # in tenths of a degree in the USHCN database. Let's go ahead and
        # correct that factor; it turns out that if you don't, the correlation
        # doesn't work correctly. Note that computing anomalies is a linear
        # operation, so it doesn't matter for the math so far that we've used
        # tenths of a degree instead of whole degrees.
        cand = series[coop_id1]
        cand_data = [val*.1 for val in cand.monthly_series]
        neighb = series[coop_id2]
        neighb_data = [val*.1 for val in neighb.monthly_series]
        
        # We SHOULD have read the same years of data, and have equal lengths
        # of data series.
        assert cand.years == neighb.years
        assert len(cand.series) == len(neighb.series)
        
        # What is the missing value placeholder? Correct for being in tenths
        # of a degree.
        MISS = cand.MISSING_VAL*.1
        
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
        r = compute_corr(cand_dif, neighb_dif, MISS)
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
    good_corrs = [coop_id2 for (coop_id2, r) in sort_corrs if r > params.corrlim]
           
    nmonths = (params.endyr-params.begyr)*12
    ksum = [0]*nmonths
    jsum = [0]*nmonths
    lowtoo = [0]*nmonths
    kstns = 0
    # Determine ksum[nmonths], the number of neighbor data available to use 
    # in homogenizing data for this station at each month
    for imo in xrange(nmonths):
        if cand_data[imo] != MISS:            
            for (k, coop_id2) in zip(xrange(len(good_corrs)), good_corrs):
                neighb = series[coop_id2]
                neighb_data = neighb.monthly_series
                if neighb_data[imo] != neighb.MISSING_VAL: 
                    ksum[imo] = ksum[imo]+1
            kstns = k
            jsum[imo] = ksum[imo]*1
                        
            if ksum[imo] < params.minpair:
                print " Total less than minpair: ",coop_id1,1900+(imo/12),1+(imo%12)
                lowtoo[imo] = 1
                
    # If we have more neighbors than necessary, then let's see if we can adjust
    # the numbers somewhat to bolster the amount of data in low-info periods,
    # being careful not too delete other good data.
    jstns = kstns*1
    if kstns > params.numcorr-1:
        
        good_corrs.reverse()
        for (k, coop_id2) in zip(xrange(len(good_corrs)), good_corrs):
            
            iremove = 1
            npair = 0
            neighb = series[coop_id2]
            neighb_data = neighb.monthly_series
            
            imonths = xrange(nmonths)
            CMISS, NMISS = MISS, neighb.MISSING_VAL
            
            iter_head = zip(imonths, cand_data, neighb_data)
            for (imo, c, n) in [(imo, c, n) for (imo, c, n) in iter_head]:
                if (c != CMISS) and (n != NMISS):
                    npair = npair+1
                    if ksum[imo] <= params.minpair:
                        print " Cannot remove:", coop_id1,'-',coop_id2,1900+(imo/12),1+(imo%12),ksum[imo],lowtoo[imo]
                        iremove = 0
                        break
                
            if iremove == 1:
                if kstns >= params.numcorr-1:
                    print " Remove:",coop_id1,"-",coop_id2,npair,corr_dict[coop_id2]
                    kstns = kstns-1
                    for imo in xrange(nmonths):
                        if cand_data[imo]!=CMISS and neighb_data[imo]!=NMISS:
                            ksum[imo] = ksum[imo]-1
        
    for imo in xrange(nmonths):
        if jsum[imo]>0:
            print "Original-Final:",coop_id1,1900+(imo/12),1+(imo%12),jsum[imo],ksum[imo]
    
    print "Original-Final Number stns:",coop_id1,jstns,kstns
                    
    all_corrs[coop_id1] = dict(corr=corr_dict)

# Write a correlation output file. The actual correlations between stations
# matches *perfectly* those computed with the ushcn_corr_2004.v3 code. However,
# there are still issues with Fortran array pointers since I don't use any here
# and I don't bother to pad the output file with bogus stations and correlations
# if there are less than we hoped to find.
print "...Assembling neighborhood correlation file\n"
corr_out = open("corr_out", 'wb')
for sta_id in test_stations:
    correlations = all_corrs[sta_id]['corr']
    sorted_neighbors = sorted(correlations.iteritems(),
                              key=itemgetter(1), reverse=True)[:params.numcorr-1]
    ids, corrs = zip(*sorted_neighbors)
    id_str = ("%6s " % sta_id)+"".join(("%6s " % id for id in ids))+"\n"
    ptr_str = ("{0: >6d} ".format(test_stations.index(sta_id)+1))+"".join(("{0: >6d} ".format(test_stations.index(id)+1) for id in ids))+"\n"
    corr_str = ("  1.00 ")+"".join(("{0: >6.2f} ".format(c) for c in corrs))+"\n"
    
    corr_out.writelines((id_str, ptr_str, corr_str))
    
    print sta_id, len(all_neighbors[sta_id])+1, len(corrs)+1
corr_out.close()

##########################################################################

# This might be something to make into a class in ushcn_data.py. As I 
# continue the workflow, I'll see what common functions I need and consider
# making a Network class.
network = dict(stations=stations, series=series, neighbors=all_neighbors,
               corrs=all_corrs)
