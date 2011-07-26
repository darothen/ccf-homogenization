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
from util import compute_mean

from mw2009.preprocess import find_neighborhood, neighborhood_strings
from mw2009.preprocess import find_correlations

from mw2009.splitmerge import diff, standardize, lrt_lookup, snht
from mw2009.splitmerge import splitmerge

from mw2009.chgptmodels import bayes, kth_line, t_test, lookup_critical

import parameters

# A Parameters object which contains the configuration for the homogenization
# analysis - specifically, the station list, settings for determining close
# neighbors, and settings for finding optimally correlated neighbors.
default_params = dict(nstns=21, 
                      mindist=200.0, 
                      distinc=200.0, 
                      numsrt=21,
                      numcorr=21, 
                      begyr=1900,
                      endyr=2001, 
                      data_src="raw", 
                      variable="avg", 
                      corrlim=0.1,
                      minpair=14,
                      benchmark=True,
                      project="benchmark")
params = parameters.default_parameters(**default_params)

pprint.pprint(params)

# Read in station data (download if necessary)
all_series, all_stations = ushcn_io.get_ushcn_data(params)

if not hasattr(params, 'stations'):
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
dist_out = open("dist_out", 'wb')

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
#for s in series.itervalues():
#    data = s.series
#    anomalies = compute_monthly_anomalies(data, -9999)
#    s.set_series(anomalies, s.years)

n.update_series = series
        
print "Determining correlated neighbors"

if os.path.exists("corr_out"):
    print "...great, I'm gonna read it from disk...",
    
    corr_file = open("data/case7.corr.ann")
    
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
    corr_out = open("corr_out", 'wb')
    for sta_id in n.stations.keys():
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
                    #stations=params.stations,
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
if os.path.exists("pair_results_benchmark"):
    print "Found pair_results on disk"
    pair_results = pickle.load(open("pair_results_benchmark", "r"))
else:
    print "Entering splitmerge to find changepoints"
    pair_results = splitmerge(n, **hom_params)
#pair_results = splitmerge(n, **hom_params)

### Try to print out the series in a table 
import numpy as np
#all_data = np.zeros([hom_params.nmo, hom_params.nstns])
station_list = ['83000%1s' % s for s in '0123456789abcdefghijk']
ids = station_list
all_data = [n.raw_series[id].monthly_series for id in ids]
all_data = np.array(all_data)
#sys.exit()
## Attribute changepoint hits to pair arrays
hits = np.zeros_like(all_data)
hits_neighbors = [[list() for i in xrange(hom_params.nmo)] for i in xrange(len(ids))]
for pair in pair_results:
    id1, id2 = pair.split("-")
    id1_ind, id2_ind = ids.index(id1), ids.index(id2)
    
    result = pair_results[pair]
    for (bp, result) in result.iteritems():
        if result['iqtype'] >= 3:
            hits[id1_ind, bp] += 1
            hits[id2_ind, bp] += 1
            
            hits_neighbors[id1_ind][bp].append(id2)
            hits_neighbors[id2_ind][bp].append(id1)
#sys.exit()
## "Unconfounding" process
#while hits.max() > 1:
#for i in [1, 2, 3, 4]:
################################################################################
#    for date_index in range(hits.shape[1]):
#        print date_index
#        ## The station/date with the highest overall changepoint count is identified
#        #max_flat_index = hits.argmax()
#        #max_indices = np.unravel_index(max_flat_index, hits.shape)
#        #assert hits[max_indices] == hits.max()
#        
#        
#        hits_this_date = hits[:,date_index]
#        #print hits_this_date.max(), hits_this_date
#        while hits_this_date.max() > 1:
#            station_index = hits_this_date.argmax()
#            #print station_index, date_index, hits_this_date[station_index]     
#            #print max_indices, hits.max()
#            
#            ## This station is the "culprit" - the cause of the breakpoints at this date
#            #station_index, date_index = max_indices
#            culprit = station_list[station_index]
#            
#            effected_neighbors = hits_neighbors[station_index][date_index]
#            ## The count on this date is decremented for all the culprit's neighbors
#            for neighbor_id in effected_neighbors:
#                neighbor_index = station_list.index(neighbor_id)
#                
#                hits[neighbor_index, date_index] -= 1
#                hits[station_index, date_index] -= 1
#    
#                hits_neighbors[neighbor_index][date_index].remove(culprit)
#                hits_neighbors[station_index][date_index].remove(neighbor_id)
#            #print hits_this_date.max(), hits_this_date
new_hits = np.zeros_like(hits)
while hits.max() > 1:
    max_in_hits = hits.max()
    max_flat_index = hits.argmax()
    max_indices = np.unravel_index(max_flat_index, hits.shape)
    
    new_hits[max_indices] = max_in_hits
    hits[max_indices] = 1    

if np.any(hits < 0): raise ValueError("hits < 0")

## Print header - 
head1 = "             |"+"|".join([i[:3] for i in ids])+"|"
head2 = "             |"+"|".join([i[3:] for i in ids])+"|"
print head1
print head2

## Print monthly series basic
from util import imo2iym
#con2str = lambda val, miss=-9999: "---" if val != miss else "-x-"
def con2str(data, missing_val=-9999):
    val, hits = data
    if val == missing_val:
        return "-x-"
    elif hits > 0:
        return "%3d" % hits
    else:
        return "---"
     
for imo in xrange(hom_params.nmo):
#for imo in xrange(10):
    year, month = imo2iym(imo)
    base_str = "%4d %2d %4d |" % (year, month, imo)
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
