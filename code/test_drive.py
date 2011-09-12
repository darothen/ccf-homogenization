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

from util import compute_monthly_anomalies
from mw2009.preprocess import find_neighborhood, neighborhood_strings
from mw2009.preprocess import find_correlations
from mw2009.preprocess import preprocess
from mw2009.splitmerge import splitmerge
from mw2009.estamt import estamt, apply_adjustments

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
                      dist_file="dist_out",
                      output_dir="data/")
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
## PRE PROCESS THE NETWORK
preprocess(n, **params)

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
                    output_dir=params.output_dir,
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

read = False
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
## COMPUTE FINAL ADJUSTMENTS FOR CHANGEPOINTS AND APPLY THEM TO THE DATA

estamt(n, **hom_params)

################################################################################
## OUTPUT SERIES
# Revert back to original data, just to be safe
station_list = n.stations.keys()
for id in station_list:
    station_series = n.raw_series[id]
    station_series.set_series(series_copy[id].series, station_series.years) 
apply_adjustments(n, **hom_params)
ushcn_io.write_network(n, **hom_params)
