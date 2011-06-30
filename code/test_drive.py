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

import os
import copy

# ccf-homogenization imports
import ushcn_io
from ushcn_data import Network
from util import compute_monthly_anomalies, scale_series
from util import diff, standardize, get_valid_data
from mw2009.preprocess import find_neighborhood, neighborhood_strings
from mw2009.preprocess import find_correlations
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
                      minpair=14,
                      project="test network")
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
for s in series.itervalues():
    data = s.series
    anomalies = compute_monthly_anomalies(data, -9999)
    s.set_series(anomalies, s.years)

n.update_series = series
        
print "Determining correlated neighbors"

if os.path.exists("corr_out"):
    print "...great, I'm gonna read it from disk...",
    
    corr_file = open("corr_out")
    
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
## BEGIN SPLITMERGE EXPERIMENTS ##

hom_params = dict(nstns=params.nstns,
                    numsrt=params.numsrt,
                    numcorr=params.numcorr,
                    begyr=params.begyr,
                    endyr=params.endyr,
                    data_src=params.data_src,
                    variable=params.variable,
                    stations=params.stations,
                    project=params.project,
                    ######################
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

id1 = "215615"
id2 = "215887"

minann = 5
begyr, endyr = hom_params.begyr, hom_params.endyr
numyr = endyr-begyr

def imo2iym(imo, begyr=1900):
    y = begyr+(imo/12)
    m = 1+(imo%12)
    return y, m

for s in series_copy.itervalues():
    data = s.series
    scaled = scale_series(data, 0.1, s.MISSING_VAL)
    anomalies = compute_monthly_anomalies(scaled, s.MISSING_VAL)
    s.set_series(anomalies, s.years)

## FIND FIRST MONTH WHERE BOTH HAVE DATA
station1 = n.stations[id1]
series1 = series_copy[id1]

data1 = series1.monthly_series
        
station2 = n.stations[id2]
series2 = series_copy[id2]

data2 = series2.monthly_series

###        
diff_data = diff(data1, data2)

MISS = series1.MISSING_VAL

## First pass through the data to find where the major "good" segment
## is. Also, will create the paired-difference series
nmo = hom_params.nmo
first = 0
first_set = False
last = 0
for (i, d1, d2) in zip(xrange(nmo), data1, data2):
    if d1!=MISS and d2!=MISS:
        if first < 12:
            first = i
            first_set = True
        last = i
inhoms = [first, last]
        
y1, m1 = imo2iym(first)
y2, m2 = imo2iym(last)

#####################################################################
## Loop over the segments we have just found and apply the 
## semihierarchical splitting algorithm and a 5% significance level
## to find where there is a split.

# first option, iopt = 1 = iTstat

segment = diff_data[first:last+1]

#z = standardize(diff_data, first, last, MISS)
z = standardize(segment, 0, numyr, MISS)

# subroutine snits - this is really just an application of the 
# likelihood ratio test. For mechanics, see Alexandersson and Moberg 1997,
# Int'l Jrnl of Climatology (pp 25-34)
#ts = [MISS for d in diff_data[first:last+1]]
ts = [MISS for d in segment]

# TODO : refactor this into util.py. Getting heavy there; might consider
# a package just for statistical tests in mw2009. 
print "Begin SNIT...",
'''
iCount = 1
for i in range(len(z)-1):
    if z[i] != MISS:
        
        zmn1, zmn2 = 0.0, 0.0
        n1, n2 = 0.0, 0.0
        for j in range(i+1):
            if z[j] != MISS:
                zmn1 = zmn1 + z[j]
                n1 = n1 + 1
                
        if n1 != 0:
            zmn1 = zmn1/n1
        else: 
            break
        
        for j in range(i+1, len(z)):
            if z[j] != MISS:
                zmn2 = zmn2 + z[j]
                n2 = n2 + 1
        
        if n2 != 0:
            zmn2 = zmn2/n2
        else:
            break
        
        ts[i] = n1*(zmn1**2) + n2*(zmn2**2)
        iCount = iCount + 1
print "done!"
'''

### NOTE -
# This is *really* counter-intuitive, and probably a bug in the original PHA
# code. Truncating the right tail of the  

iCount = 1
mcnt = len(get_valid_data(z))
#for pivot in range(len(z)-1):
for pivot in range(mcnt-1):
    if z[pivot] != MISS:
        
        left_series = get_valid_data(z[:pivot+1])
        #right_series = get_valid_data(z[pivot+1:])
        right_series = get_valid_data(z[pivot+1:mcnt])
        
        zmleft = sum(left_series)
        nleft = len(left_series)
        if nleft != 0:
            zmleft = zmleft/nleft
        else: 
            break
        
        zmright = sum(right_series)
        nright = len(right_series)
        if nright != 0:
            zmright = zmright/nright
        else:
            break
        
        ts[pivot] = nleft*(zmleft**2) + nright*(zmright**2)
        iCount = iCount+1
print "...done again!"
        
# Okay, we have likelihood ratios at each segment. Now, find the maximum
# one and its index
iPeak = 0
rPeak = 0.0
clip_ts = ts[2:-2] # We clip the beginning and end
for (ind, ts_val) in zip(xrange(len(clip_ts)), clip_ts):
    if ts_val > rPeak:
        iPeak = ind
        rPeak = ts_val
# Now we find the critical value for this data, and check our max
# likelihood ratio against it
def t_lookup(num_vals):
    '''alpha = 0.5
    
    Alexandersson and Moberg 1997,
    Int'l Jrnl of Climatology (pp 25-34)
    Appendix 2, Table AI
    '''
    sig95 = [4.54,5.70,6.95,7.65,8.10,8.45,8.65,8.80,8.95,9.05,
             9.15,9.35,9.55,9.70]
    nsig = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250]
    
    if num_vals < nsig[0]: 
        return 99999.0 # Too few values!
    elif num_vals >= nsig[-1]:
        return sig95[-1] # This is a good approximation as n->infinity
    else:
        # Select the two values from nsig which bound num_vals
        bounds = zip(nsig[:-1], nsig[1:])
        bi = 0
        l, r = bounds[bi]
        while not (l <= num_vals <= r): 
            bi = bi+1
            l, r = bounds[bi]
        left_ind, right_ind = nsig.index(l), nsig.index(r)
        
        # Estimate the critical value using linear interpolation between
        # the two lookup values we found
        bound_left, bound_right = l, r
        t_left, t_right = sig95[left_ind], sig95[right_ind]
        
        interp = num_vals - bound_left
        slope = (t_right - t_left) / (bound_right - bound_left)
        
        t_val = t_left + slope*interp
        return t_val
    
crit_val = t_lookup(iCount)
# Test the peak stat against the critical value
curstat = rPeak
new_end_mo = first + iPeak + 2 # shift by 2 aligns clip_ts with ts,
                               # shift by first aligns it back to diff_data
ynew, mnew = imo2iym(new_end_mo)
# Fragment First if either homog or inhomog
print "%6s-%6s MD       FIRST series %4d %2d to %4d %2d | at %4d %2d ts: %4.2f limit >: %3.2f" % (id1,id2,y1,m1,y2,m2,ynew,mnew,curstat,crit_val)
print z[0], ts[0]
