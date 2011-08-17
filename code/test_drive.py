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

from mw2009.estamt import estamt

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
default_params = dict(nstns=52, 
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
                    minlenshf=24,
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
                    qscale=1.46 # scaling factor for clipping outliers
                )
hom_params = parameters.Parameters(**hom_params)

import pickle
pairs = list(pickle.load(open("fortran_pairs", "r")))
hom_params['pairs'] = pairs

read = True
if os.path.exists("pair_results") and read:
    print "Found pair_results on disk"
    pair_results = pickle.load(open("pair_results", "r"))
else:
    print "Entering splitmerge to find changepoints"    
    pair_results = splitmerge(n, **hom_params)
n.pair_results = pair_results

from mw2009.confirmfilt import filter1, filter2, filter3
filter1(n, **hom_params)
filter2(n, **hom_params)
filter3(n, **hom_params)

#################################################################################
## COMPUTE FINAL ADJUSTMENTS FOR CHANGEPOINTS

estamt(n, **hom_params)

################################################################################
## OUTPUT SERIES
print "----------- Output Adjustments -----------"
out_dir = "data/WMs.52d"

flags = ' ABCDEFGHIJKLMNOPQRSTUVWXYZ'
station_list = n.stations.keys()

#station_list = ['215887', ]
for id in station_list:    
    station_series = n.raw_series[id]
    #station_data = n.raw_series[id].monthly_series[:]
    station_data = series_copy[id].monthly_series[:]
    station_index = station_list.index(id)
    miss = n.raw_series[id].MISSING_VAL

    station_changepoints = n.raw_series[id].changepoints
    cps = sorted(station_changepoints.keys())
    nchg = len(cps)
    first, last = cps[0], cps[-1]
    
    adjusted_data = [miss]*len(station_data)
    adjtemp = [miss]*len(station_data)
    contemp = [miss]*len(station_data)
    outtemp = [miss]*len(station_data)
    adjflag = [' ']*len(station_data)
    
    sumchg = 0.0
    sumcon = 0.0
    jadj = 0
    aflg = ' '
    
    segs_to_adjust = zip(cps[-2::-1], cps[::-1])
    for left, right in segs_to_adjust:
        jadj += 1
                
        cp_index = cps.index(left)
        print "imo: ",station_index,cp_index,(left+1),right
        for imo in range(right, left, -1):
            if station_data[imo] != miss:
                totchg = sumchg
                adjtemp[imo] = totchg
                contemp[imo] = sumcon
                outtemp[imo] = station_data[imo]*.1 - totchg
                adjflag[imo] = aflg
                
        print ("Adj write: %s %5d %5d %5d %5d %5d %7.2f %7.2f %s" % 
               (id,nchg,left,right,cp_index,jadj,sumchg,sumcon,aflg) )
        
        adj = station_changepoints[left]['ahigh']
        std = station_changepoints[left]['astd']
        
        aflg = flags[jadj]
        sumchg = totchg+adj
        sumcon = sqrt(sumcon**2 + std**2)
        
    ## call writsta(itarg, nstn(itarg), outtemp, adjtemp, contemp, adjflag, otag, idunit)
    for i in range(len(outtemp)):
        o = outtemp[i]
        if o != miss:
            outtemp[i] = int(o*10.)
            
    yearly_outtemp = []
    for y in range(hom_params.numyr):
        yearly_outtemp.append(outtemp[y*12:(y+1)*12] + [miss])
    station_series.set_series(yearly_outtemp, station_series.years)
        
    out_filename = "%s_%s.WMs.52d" % (id, hom_params.variable)
    output_file = open(os.path.join(out_dir, out_filename), 'w')
    print " Writing: CoopOutDir:",out_filename
    output_lines_to_write = []
    for (year_data, year) in zip(station_series.series, station_series.years):
        write_data = "".join(["%6d" % val for val in year_data])
        outstr = "%6s1%4d%s\n" % (id,year,write_data)
        output_lines_to_write.append(outstr)
    output_file.writelines(output_lines_to_write)
    output_file.close()

sys.exit()