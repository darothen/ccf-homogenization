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

import pickle

from math import sqrt

import os
import sys
import copy

# ccf-homogenization imports
import ushcn_io
from ushcn_data import Network

from util import compute_monthly_anomalies, scale_series, get_valid_data
from util import compute_mean, imo2iym, compute_monthly_avg_std

from mw2009.preprocess import find_neighborhood, neighborhood_strings
from mw2009.preprocess import find_correlations

from mw2009.splitmerge import diff, standardize, lrt_lookup, snht
from mw2009.splitmerge import splitmerge

from mw2009.chgptmodels import bayes, kth_line, t_test, lookup_critical

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
                      numcorr=40, 
                      begyr=1900,
                      endyr=2010, 
                      data_src="raw", 
                      variable="max", 
                      stations=test_stations,
                      corrlim=0.1,
                      minpair=14,
                      project="test network",
                      corr_file="corr_out_40",
                      dist_file="dist_out")
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

##
series_copy = copy.deepcopy(series)
##

n = Network(stations, series, name=params.project)
print n

##########################################################################
print "Analyzing geographic network neighborhoods"

all_neighbors = dict()
stations_list = n.stations.values()
for station in stations_list:
    
    print station.coop_id
    print "...computing neighbor distances"
    
    neighbors = find_neighborhood(station, stations_list, **params)
    all_neighbors[station.coop_id] = neighbors
    
n.neighborhoods = all_neighbors
 
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
dist_out = open(params.dist_file, 'wb')

for station in stations_list:
    print "   ", station
    neighbors = all_neighbors[station.coop_id]
    
    out_strings = neighborhood_strings(station, neighbors, stations_list)
    dist_out.writelines(out_strings)

dist_out.close()
    
##########################################################################

# Go through all the data we have, and replace the read-in values with
# monthly anomalies. Then, flatten the data into a list with all the data
# and length (endyr-begyr)*12
for s in series.itervalues():
    data = s.series
    anomalies = compute_monthly_anomalies(data, -9999)
    s.set_series(anomalies, s.years)

n.update_series = series
        
print "Determining correlated neighbors"

if os.path.exists(params.corr_file):
    print "...great, I'm gonna read it from disk...",
    
    corr_file = open(params.corr_file)
    
    all_lines = corr_file.readlines()
    station_lines = all_lines[::2]
    corr_lines = all_lines[1::2]
    
    all_corrs = dict()
    for (sta_line, corr_line) in zip(station_lines, corr_lines):
        stations = sta_line.strip().split()
        this, others = stations[0], stations[1:]
        corrs = map(float, corr_line.strip().split()[1:])
        
        corr_dict = dict()
        for (id, corr) in zip(others, corrs):
            if not id == '000000':
                corr_dict[id] = corr
        
        all_corrs[this] = corr_dict
        
    print " that was fast!"
    
else:
    
    print "...shit, need to do all the compuations"

    all_corrs = dict()
    
    for cand_series in n.raw_series.itervalues():
        
        coop_id1 = cand_series.coop_id
        corr_dict = find_correlations(cand_series, n.raw_series,
                                      n.neighborhoods[coop_id1], **params)
        
        all_corrs[coop_id1] = dict(corr=corr_dict)
    
    # Write a correlation output file. The actual correlations between stations
    # matches *perfectly* those computed with the ushcn_corr_2004.v3 code. However,
    # there are still issues with Fortran array pointers since I don't use any here
    # and I don't bother to pad the output file with bogus stations and correlations
    # if there are less than we hoped to find.
    print "...Assembling neighborhood correlation file\n"
    corr_out = open(params.corr_file, 'wb')
    for sta_id in test_stations:
        correlations = all_corrs[sta_id]['corr']
        sorted_neighbors = sorted(correlations.iteritems(),
                                  key=itemgetter(1), reverse=True)[:params.numcorr-1]
        # PAD PAD PAD
        ## Just pad the output with '000000' stations, r = 0.0 to make it look
        ## like the normal output.
        while len(sorted_neighbors) < params.numcorr-1:
            sorted_neighbors.append(('000000', 0.0))
        ids, corrs = zip(*sorted_neighbors)
        id_str = ("%6s " % sta_id)+"".join(("%6s " % id for id in ids))+"\n"
        #ptr_str = ("{0: >6d} ".format(test_stations.index(sta_id)+1))+"".join(("{0: >6d} ".format(test_stations.index(id)+1) for id in ids))+"\n"
        corr_str = ("  1.00 ")+"".join(("{0: >6.2f} ".format(c) for c in corrs))+"\n"
        
        #corr_out.writelines((id_str, ptr_str, corr_str))
        corr_out.writelines((id_str, corr_str))
        
    corr_out.close()

n.correlations = all_corrs

##########################################################################
## BEGIN SPLITMERGE EXPERIMENTS ##if os.path.exists("corr_out"):

hom_params = dict(nstns=params.nstns,
                    numsrt=params.numsrt,
                    numcorr=params.numcorr,
                    beg_year=params.begyr,
                    end_year=params.endyr,
                    data_src=params.data_src,
                    variable=params.variable,
                    stations=params.stations,
                    project=params.project,
                    ######################if os.path.exists("corr_out"):
                    numyr=params.endyr-params.begyr,
                    nmo=(params.endyr-params.begyr)*12,
                    minser=5, # min number of ind. months in a raw series that can be tested
                    minann=5, # min number of years for a given month
                    slpmthmax=0.0225, # max slope (degrees/month) threshold for sloped models
                    minlen=18, # min number of (months) for an estimate window
                    mincomp=60, #
                    compt=0.8, # if segment has more than minlen but less than mincomp,
                               # then compt is % completeness needed for estimate
                    minsta=2, # min number of station pairs for estimate eval
                    minhits=2, # min num of sta-neigh series at the model decision phase
                                # required to use interquartile estimate, else this is a 
                                # linear (no chgpt) model with amp := 0.0
                    homog=1,   # --
                    indeter=2, # test statistics decisions by all of the techniques
                    inhomog=3, # --
                    ninh=80, # max number of breaks in the series
                    inhnet=440, # max number of breaks in a network
                    eps=1e-6, # a very small number
                    stepthres=0.0, # a temperature step limit at which these models might work
                )
hom_params = parameters.Parameters(**hom_params)

import pickle
pairs = list(pickle.load(open("fortran_pairs", "r")))
#pairs = [('314684', '315177'),
#         ('086997', '314684'), 
#         ('307633', '314684'),
#         ('124837', '324418'),
#         ('164407', '314055')]
#pairs = [('111280', '124837'), ]

read = True
if os.path.exists("pair_results") and read:
    print "Found pair_results on disk"
    pair_results = pickle.load(open("pair_results", "r"))
else:
    print "Entering splitmerge to find changepoints"    
    pair_results = splitmerge(n, **hom_params)
    
## TEMP TEMP TEMP - remove slr hits from pair_results:
for pair_str in pair_results.keys():
    
    for bp in pair_results[pair_str].keys():
        if bp == 'del': 
            continue
        
        delete_bps = []
        if 'SLR' in pair_results[pair_str][bp]['cmodel']:
            delete_bps.append(bp)
        for bp in delete_bps:
            del pair_results[pair_str][bp]

### Try to print out the series in a table 
import numpy as np
#all_data = np.zeros([hom_params.nmo, hom_params.nstns])
station_list = n.stations.keys()
ids = station_list
all_data = [n.raw_series[id].monthly_series for id in ids]
all_data = np.array(all_data)

## Attribute changepoint hits to pair arrays
hits = np.zeros_like(all_data)

hits_neighbors = [[list() for i in xrange(hom_params.nmo)] for i in xrange(len(ids))]
ntests = np.zeros_like(all_data)
#for pair in pairs:
#    id1, id2 = pair
#    ntests[ids.index(id1), :] += 1
#    ntests[ids.index(id2), :] += 1

for pair in pair_results:
    id1, id2 = pair.split("-")
    id1_ind, id2_ind = ids.index(id1), ids.index(id2)
    
    ntests[id1_ind, :] += 1
    ntests[id2_ind, :] += 1
        
    result = pair_results[pair]
    for (bp, result) in result.iteritems():
        if bp == 'del':
            #continue
            print pair, pair_results[pair]['del']
            for bad_month in pair_results[pair]['del']:
                hits[id1_ind, bad_month] += 1
                hits[id2_ind, bad_month] += 1
            
        else:
            #if 'TPR' in result['cmodel']:
            hits[id1_ind, bp] += 1
            hits[id2_ind, bp] += 1
       
            hits_neighbors[id1_ind][bp].append(id2)
            hits_neighbors[id2_ind][bp].append(id1)
            
################################################################################
## PRE FILTER 1

################################################################################
## FILTER 1 TEST
new_hits = np.zeros_like(hits)
amps = np.zeros_like(hits)

hit_tuples = []

for imo in range(hits.shape[1]):
    hits_month = hits[:,imo]
    while hits_month.max() > 1:
        
        max_in_hits = hits_month.max()
        station_index = hits_month.argmax()
    
        station_id = station_list[station_index]
        iy, im = imo2iym(imo)
        
        #print max_in_hits, station_index, station_id
        
        ## Find entries in pair_results with this station and a changepoint on this
        ## date
        pr_keys = [key for key in pair_results if (station_id in key and
                                                   imo in pair_results[key])]
        #print pr_keys
        
        if not pr_keys: 
            break
        offset_sum, offset_z_sum, count = 0.0, 0.0, 0.0
        for key in pr_keys:
            bp_summary = pair_results[key][imo]

            id1, id2 = key.split("-")
            
            offset, offset_z = bp_summary['offset'], bp_summary['offset_z']
            if station_id == id2: 
                offset = offset*-1.0
            
            print ntests[station_index, imo],"->",
            ntests[station_index, imo] -= 1
            print ntests[station_index, imo]
            
                
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
          
        #ntests[station_index, imo] -= 1
        
        hit_tuples.append( (station_index, imo) )
        
# set ntests
#for station_index, month in hit_tuples:
#    ntests[station_index,        :month] -= 1
#    ntests[station_index, month+1:     ] -= 1
        
if np.any(hits < 0): raise ValueError("hits < 0")

################################################################################
## FILTER 2 TEST
## won't actually do anything but print things to console for now
print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
print "FILTER 2 - SHF"
print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
old_station_data = n.raw_series['111280'].series[:2]
(num_stations, num_months) = amps.shape
for station_index in range(num_stations):
    
    station_id = station_list[station_index]
    data = n.raw_series[station_id].series
    scale_series(data, 0.1, s.MISSING_VAL)
    stdk = compute_monthly_avg_std(data)
    
    for month in range(num_months):
        
        if new_hits[station_index, month] > 0:
            ahigh = amps[station_index, month]
            jsum = new_hits[station_index, month]
            
            astd = abs(ahigh)
            
            iy, im = imo2iym(month)
            
            ## This are PRE-DEFINED in inhomog.parm.system.mthly.incl. They 
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

################################################################################
## FILTER 3 TEST           
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

inconfirm = 2
nmrgyr = -2

final_hits = np.zeros_like(hits)
ifound = 0

print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
print "FILTER 3 - COLLAPSE CHANGEPOINT SAMPLES"
print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
(num_stations, num_months) = amps.shape
for station_index in range(num_stations):
    
    station_id = station_list[station_index]
    data = n.raw_series[station_id].series
    miss = n.raw_series[station_id].MISSING_VAL
    data_monthly = n.raw_series[station_id].monthly_series
    stdk = compute_monthly_avg_std(data)
    
    ## setup temp arrays 
    khits = np.zeros(hom_params.nmo)
    ktests = np.zeros(hom_params.nmo)
    akhigh = np.zeros(hom_params.nmo)
    
    # iterate until there are no more high points
    istop = False
    while not istop:
        
        ihighit, ahigh = 0.0, 0.0
        
        # find the highest count - the most number of hits at a possible breakpoint
        for month in range(num_months):
            
            isum, asum = 0.0, 0.0
            if new_hits[station_index, month] >= inconfirm:
                jhit = new_hits[station_index, month]
                isum += jhit
                asum += amps[station_index, month]*jhit
            
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
                if right_month == hom_params.nmo:
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
                    new_hits[station_index,ihighmo] = 0
                    amps[station_index,ihighmo] = 0.0
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
                    new_hits[station_index,ihighmo] = 0
                    amps[station_index,ihighmo] = 0.0
                    absorbed = True
                    break
            
                radius += 1
                    
            # if no hits found, setup new hit
            if not absorbed:
                khits[ihighmo] = ihighit
                ktests[ihighmo] = ntests[station_index,ihighmo]
                akhigh[ihighmo] = ahigh*ihighit
                print "New CHG hit: ",station_index,ihighmo,khits[ihighmo],ktests[ihighmo],akhigh[ihighmo]/khits[ihighmo]
        
                new_hits[station_index, ihighmo] = 0
                amps[station_index, ihighmo] = 0.0
            
            #raw_input("pause")
        else:
            istop = True
            
    print "----------------------------------------------"
    # examine interim khits array for station's filtered changepoints
    for month in range(hom_params.nmo):
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
    
print "-------------------------------------------------"
print "Undoc filter: ",ifound
sys.exit()

## Print header - 
head1 = "             |"+"|".join([i[:3] for i in ids])+"|"
head2 = "             |"+"|".join([i[3:] for i in ids])+"|"
print head1
print head2

## Print monthly series basic
#con2str = lambda val, miss=-9999: "---" if val != miss else "-x-"
def con2str(data, missing_val=-9999):
    val, hits = data
    
    if hits > 0:
        return "%3d" % hits
    #elif val == missing_val:
    #    return "-X-"
    else:
        return "---"
     
for imo in xrange(hom_params.nmo):
#for imo in xrange(10):
    year, month = imo2iym(imo)
    base_str = "%4d %2d %4d |" % (year, month, imo-11)
    month_strs = "|".join(map(con2str, zip(all_data[:,imo],
                                           new_hits[:,imo])))+"|"
    print_month_strs = False
    for i in range(10):
        if str(i) in month_strs: 
            print_month_strs = True
            break
    
    if print_month_strs:
        print base_str+month_strs
    
# Test plot
#pair_str = pair_results.keys()[-1]
#id1, id2 = pair_str.split("-")
#bp_inds = pair_results[pair_str].keys()
#bps = pair_results[pair_str]
#bp_colors = [bps[bp]['iqtype'] for bp in bp_inds]
#bp_cs = []
#for c in bp_colors:
#    if c > 3:
#        bp_cs.append('r')
#    else:
#        bp_cs.append('y')
#
#series1, series2 = n.raw_series[id1], n.raw_series[id2]
#
#anom1 = compute_monthly_anomalies(series1.series, -9999)
#series1.set_series(anom1, series1.years)
#data1 = series1.monthly_series
#
#anom2 = compute_monthly_anomalies(series2.series, -9999)
#series2.set_series(anom2, series2.years)
#data2 = series2.monthly_series
#
#import numpy as np
#from pylab import *
#
#dm1 = np.ma.masked_equal(data1, series1.MISSING_VAL)*.1
#dm2 = np.ma.masked_equal(data2, series2.MISSING_VAL)*.1
#
#subplot(2,1,1)
#plot(dm1)
#vlines(bp_inds, ylim()[0], ylim()[1], colors=bp_cs)
#title(series1)
#
#subplot(2,1,2)
#plot(dm2)
#vlines(bp_inds, ylim()[0], ylim()[1], colors=bp_cs)
#title(series2)
#
#subplots_adjust(hspace=.3)
