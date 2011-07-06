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

from math import sqrt

import os
import copy

# ccf-homogenization imports
import ushcn_io
from ushcn_data import Network

from util import compute_monthly_anomalies, scale_series, get_valid_data
from util import compute_mean

from mw2009.preprocess import find_neighborhood, neighborhood_strings
from mw2009.preprocess import find_correlations

from mw2009.splitmerge import diff, standardize, lrt_lookup, snht

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

#id1 = "215887"
#id2 = "215615"
id2 = "153430"
id1 = "034572"

minann = 5
begyr, endyr = hom_params.begyr, hom_params.endyr
numyr = endyr-begyr

def imo2iym(imo, begyr=1900):
    y = begyr+(imo/12)
    m = 1+(imo%12)
    return y, m
    
def within(test_interval, cover):
    
    test_l, test_r = test_interval
    cover_l, cover_r = cover
    return (cover_l <= test_l <= test_r <= cover_r)

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

breakpoints = []
breakpoints.extend(inhoms)
homog_segs = []

iter = 0
enter_BIC = False
last_breakpoints = []
while (iter < 10) and not enter_BIC:
    
    #print map(imo2iym, breakpoints) 
    
    seg_bounds = zip(breakpoints[:-1], breakpoints[1:])
    last_breakpoints = copy.deepcopy(breakpoints)
    new_breakpoints = copy.deepcopy(breakpoints)
        
    seg_lookup = []
    
    new_homog_segs = []

    new_test_segs = []

    print "Parse segments (isplit = 1), ipass: ", iter
    print breakpoints
    print seg_bounds
    for (l, r) in seg_bounds:
                
    #####################################################################
    ## Loop over the segments we have just found and apply the 
    ## semihierarchical splitting algorithm and a 5% significance level
    ## to find where there is a split.
        y1, m1 = imo2iym(l)
        y2, m2 = imo2iym(r)

        adjust = int(seg_bounds.index((l,r)) > 0)
        segment = diff_data[l+adjust:r+1]
        numyr = r-l
        
        # short circuit and skip if it's too short
        if (r-l) <= 5:
            print "Too short: ", imo2iym(l), imo2iym(r)
            continue
        
        # short circuit and skip this segment if we already know that it's
        # homogenous
        this_within = lambda seg: within((l, r), seg)
        within_stable_segs = map(this_within, homog_segs)
        if any(within_stable_segs):
            print "Stable segment: ", imo2iym(l), imo2iym(r)
            if l == first: 
                new_breakpoints.append(first)
            seg_lookup.append(((l, r), 'stable'))
            continue
        
        z = standardize(segment, MISS)
        
        # subroutine snits - this is really just an application of the 
        # likelihood ratio test. For mechanics, see Alexandersson and Moberg 1997,
        # Int'l Jrnl of Climatology (pp 25-34)
        ts = snht(z, MISS, standardized=True)
        z_count = len(get_valid_data(z))
                
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
        crit_val = lrt_lookup(z_count)
        # Test the peak stat against the critical value
        curstat = rPeak
        new_end_mo = l + iPeak + 2 + adjust # shift by 2 aligns clip_ts with ts,
                                            # shift by first aligns it back to diff_data
                                            # shift by adjust aligns it back with the diff_data
        ynew, mnew = imo2iym(new_end_mo)
        # Fragment First if either homog or inhomog
        if iter == 0: # first time through loop
            print "%6s-%6s MD        FIRST series %4d %2d to %4d %2d | at %4d %2d ts: %4.2f limit >: %3.2f" % (id1,id2,y1,m1,y2,m2,ynew,mnew,curstat,crit_val)
            breakpoints.append(new_end_mo)
            breakpoints = sorted(breakpoints)
            seg_status = 'first'
            
        else:
            if curstat > crit_val:
                print "%6s-%6s MD Inhomogenity for series %4d %2d to %4d %2d | at %4d %2d ts: %4.2f limit >: %3.2f %4d" % (id1,id2,y1,m1,y2,m2,ynew,mnew,curstat,crit_val,z_count)
                new_breakpoints.append(new_end_mo)
                #seg_status = 'inhomog'
                #seg_lookup.append(((l, r), 'inhomog', new_end_mo))
                        
            else:
                print "%6s-%6s MD      Homogeneous series %4d %2d to %4d %2d | at %4d %2d ts: %4.2f limit >: %3.2f %4d" % (id1,id2,y1,m1,y2,m2,ynew,mnew,curstat,crit_val,z_count)
                #seg_lookup.append(((l, r), 'homog'))
                #seg_status = 'homog'
                
                new_homog_segs.append((l, r))
                
    # update our account of which segments were homogenous. We want to remember
    # them during the next iteration. Also, condense stable segments if possible
    # (i.e., two lie next to each other, replace them with the union of those
    # two segments)
    homog_segs.extend(new_homog_segs)
    if homog_segs:
        homog_segs = sorted(homog_segs, key=itemgetter(0))
        final_homog_segs = [homog_segs[0], ] # this will be like a stack
        for seg in homog_segs[1:]:
            last_seg = final_homog_segs[-1]
            if last_seg[1] == seg[0]:
                new_seg = (last_seg[0], seg[1])
                final_homog_segs.pop()
                final_homog_segs.append(new_seg)
            else:
                final_homog_segs.append(seg)
        homog_segs = final_homog_segs

    # So we have new segments that can be generated from these new
    # breakpoints. Now, the PHA routine enters a "merge" process
    # to see whether or not to keep these newly found changepoints or throw
    # them out as false alarms. The basis of this is whether or not the 
    # test statistic we found from the likelihood ratio test is significant
    # or not. Now, we will double check them, and if it is not significant,
    # we will remove that breakpoint and reapeat the splitting process.
    new_breakpoints = sorted(new_breakpoints)
    print new_breakpoints
    seg_bounds = zip(new_breakpoints[:-2], new_breakpoints[2:])
    print seg_bounds
    
    remove_breakpoints = set()
    add_breakpoints = set()
    merged_breakpoints = set()
    if iter > 0:
        
        collapse_segs = []
        
        print "Merge segments (isplit = 0), ipass: ", iter
        
        # An annoying thing here is that we will need to find potential breakpoints 
        # between the new ones we found. Let's get that out of the way:
        #seg_bounds = merge_breakpoints
        
        #for (l, r) in seg_bounds:
        for ((l ,r), new_bp) in zip(seg_bounds, new_breakpoints[1:-1]):
            
            print imo2iym(l), imo2iym(r)
            
            # short circuit and skip this segment if we already know that it's
            # homogenous
#            this_within = lambda seg: within((l, r), seg)
#            within_stable_segs = map(this_within, homog_segs)
#            if any(within_stable_segs):
#                print "Stable segment: ", imo2iym(l), imo2iym(r)
#                if l == first: 
#                    new_breakpoints.append(first)
#                seg_lookup.append(((l, r), 'stable'))
#                continue
            
            adjust = int(seg_bounds.index((l, r)) > 0)
            segment = diff_data[l+adjust:r+1]
            
            z = standardize(segment, MISS)
            ts = snht(z, MISS, standardized=True)
            z_count = len(get_valid_data(z))
                
            iPeak = 0.0
            rPeak = 0.0
            clip_ts = ts[2:-2] # We clip the beginning and end
            for (ind, ts_val) in zip(xrange(len(clip_ts)), clip_ts):
                if ts_val > rPeak:
                    iPeak = ind
                    rPeak = ts_val
            crit_val = lrt_lookup(z_count)
            curstat = rPeak
            new_end_mo = l + iPeak + 2 + adjust
            
            iy, im = imo2iym(new_end_mo)
            
            if curstat > crit_val:
                print "%6s-%6s MD  Peak kept in merge at %4d %2d | ts: %4.2f limit >: %3.2f" % (id1,id2,iy,im,curstat,crit_val)

                merged_breakpoints.add(l)
                merged_breakpoints.add(new_bp)
                merged_breakpoints.add(r)
            else:
                print "%6s-%6s MD Compress 2 out peak at %4d %2d | ts: %4.2f limit >: %3.2f" % (id1,id2,iy,im,curstat,crit_val)
                # Crap, if there are any potential breakpoints in this segment,
                # we need to remove them because this segment is homogenous. Let's
                # remember this homogenous segment for now and come back once
                # we've found all of them.

                merged_breakpoints.update([l, r])
                remove_breakpoints.add(new_bp)

        merged_breakpoints.difference_update(remove_breakpoints)
        breakpoints = list(merged_breakpoints)
    breakpoints = sorted(breakpoints)
    
    enter_BIC = (breakpoints == last_breakpoints)
    iter = iter + 1
    
## Okay wow, we've potentially made it to the BIC stage now... !
if first not in breakpoints:
    breakpoints.insert(0, first)
ym_breakpoints = map(imo2iym, breakpoints)
print ym_breakpoints

## ENTERING MINBIC
#left, bp, right = breakpoints[6:9]
left, bp, right = breakpoints[1:4]
if left != first:
    left = left + 1
# recall that we only consider data after the first full year. we will be 
# computing regressions with the independent variable indexed from this 
# starting point, so we need to shift these indices. we also need to shift them
# by +1 if this is any segment beyond the first one, so that we don't include
# changepoints in more than one analysis.
# TOTAL_SHIFT = -12 + 1 = -11
# 
# However, this shift is only necessary while looking at the array indices that
# we generate using range(). the data should already be aligned correctly.
total_shift = -12 + 1
left_shift, bp_shift, right_shift = left+total_shift, bp+total_shift, right+total_shift
y1, m1 = imo2iym(left)
yb, mb = imo2iym(bp)
y2, m2 = imo2iym(right)
print "Entering MINBIC - %4d %2d    %4d %2d    %4d %2d" % (y1, m1, yb,
                                                           mb, y2, m2)

## Print header for BIC changepoint testing -
left_header = " QTYP     QVAL    QRSE     QPF     MU1     MU2  ALPHA1"
right_header = "  ALPHA2   MSTAT   MCRIT    MOFF KNT1 KNT2"
print (left_header+right_header)

################################################################################     

## Looks like the first test is KTHSLR1, kendall-theil method with slope calc
## We perform this test on the entire interval containing the breakpoint
(seg_x, seg_data) = range(left_shift, right_shift+1), diff_data[left:right+1] 
cmodel = "KTHSLR1"
#lsql = least_squares(seg_x, seg_data, MISS)
kthl = kth_line(seg_x, seg_data, MISS)

nobs = right-left+1
slpmed = kthl.slope
yintmed = kthl.y_int
sseredmed = kthl.sseslope
nval = kthl.nval
qoff = 0.0

qslr1, rsq1, rsq2 = bayes(nobs, sseredmed, 2)
# output string
head = "%7s %6.2f %7.2f %7.2f" % (cmodel, qslr1, rsq1, rsq2)
stats =  " %7.2f ------- %7.3f ------- ------- -------" % (yintmed, slpmed)
tail = "% 7.2f %5d ----" % (qoff, nval)
print (head+stats+tail)

################################################################################     

## Now we begin the two-phase regressions, where we use the kendall-theil fit
## on the data at either side of the breakpoint.
##
## The first regression assumes that the segments have the same slope
left_seg = range(left_shift, bp_shift+1)
left_data = diff_data[left:bp+1]
right_seg = range(bp_shift+1, right_shift+1)
right_data = diff_data[bp+1:right+1]

# kendall-theil method with 0 sloped segments
cmodel = "KTHTPR0"

kthl_left = kth_line(left_seg, left_data, MISS)
kthl_right = kth_line(right_seg, right_data, MISS)

left_y_med = kthl_left.y_med
right_y_med = kthl_right.y_med
q_off = left_y_med - right_y_med

n_left = kthl_left.nval
n_right = kthl_right.nval
n_total = n_left + n_right

stat_test = t_test(left_data, right_data, MISS)
t_val = stat_test.t_val
t_crit = lookup_critical(n_total-2, "t")

sse_sum = kthl_left.sseflat + kthl_right.sseflat
qslr1, rsq1, rsq2 = bayes(n_total, sse_sum, 3)
# output string
head = "%7s %6.2f %7.2f %7.2f" % (cmodel, qslr1, rsq1, rsq2)
stats =  " %7.2f %7.2f ------- ------- %7.2f %7.2f" % (left_y_med, right_y_med,
                                                      t_val, t_crit)
tail = "% 7.2f %5d %4d" % (q_off, n_left, n_right)
print (head+stats+tail)

################################################################################     

## The second regression tests for a step change with equal (constant) sloped
## segments

all_data = diff_data[left:right+1]
nobs = right+1-left
left_data = diff_data[left:bp+1]
right_data = diff_data[bp+1:right+1]

## 1) Compute the mean for *all* of the data
all_valid_data = get_valid_data(all_data, MISS)
all_mean = compute_mean(all_valid_data, valid=True)

## 2) use kendall-theil method with single slope
# This method is slightly different than the kth_line() method above.
# First, we get only the valid data, and we pair it with the natural ordering
# of the data, i.e. 1, 2, 3...
valid_all = all_valid_data
n_all = len(valid_all)
range_all = range(1, n_all+1)

valid_left = get_valid_data(left_data, MISS)
valid_right = get_valid_data(right_data, MISS)

n_left, n_right = len(valid_left), len(valid_right)
range_left = range(1, n_left+1)
range_right = range(n_left+1, n_all+1)

# Second, generate paired slopes for the first segment.
nslp = 0
r_temp = []
for i in range(n_left-1):
    for j in range(i, n_left):
        if range_left[j] != range_left[i]:
            nslp = nslp + 1
            r_temp.append( (valid_left[j]-valid_left[i])/
                           (range_left[j]-range_left[i]) )
# Third, generate paired slopes for the second segment.
for i in range(n_right-1):
    # BUG: MW2009 code in chgptmodels.kendallthiell, line 2229 starts the 'j'
    #     index at ibeg2+1. This corresponds to 1 here. Above in the first
    #     segment and in kth_line(), it starts the 'j' index right where 'i' 
    #     left off.
    for j in range(1, n_right):
        if range_right[j] != range_right[i]:
            nslp = nslp + 1
            r_temp.append( (valid_right[j]-valid_right[i])/
                           (range_right[j]-range_right[i]) )
            
#Fourth, find the median slope from all the ones we computed
islope = 1
if not islope:
    r_slope = 0.0
else:
    r_temp = sorted(r_temp)
    imed = (nslp - 1)/2
    if (nslp%2)==1: imed = imed+1 # offset by one to right if odd
    r_slope = r_temp[imed]
    
print "slope, ic, imet: %7.2f %5d %5d" % (r_slope, nslp, imed)

# Fifth, compute the first segment intercept, y-median - slope*x-median
imed = (n_left - 1)/2
if (n_left%2)==1: imed = imed+1
range_med = range_left[imed]
valid_left = sorted(valid_left)
data_med = valid_left[imed]
left_y_int = data_med-r_slope*range_med
print "Seg1 - Xmed, Ymed, slope, Yint: %7.2f %7.2f %7.2f %7.3f" % (range_med, 
                                                                   data_med,
                                                                   r_slope,
                                                                   left_y_int)
# BUG: Again in chgptmodel.kendalltheill(), there is a bug on line 2339. Starting
#     here, we over-write the medians we found in both lists, and use the second
#     segment for all our computations! I reproduce that behavior ehre by using the
#     generic range_med and data_med values for rXmed and rYmed. Should we not
#     care about different medians for each segment?
# Sixth, compute the second segment intercept
imed = (n_right - 2)/2
if (n_right%2)==1: imed = imed + 1
range_med = range_right[imed]
valid_right = sorted(valid_right)
data_med = valid_right[imed]
right_y_int = data_med-r_slope*range_med
print "Seg2 - Xmed, Ymed, slope, Yint: %7.2f %7.2f %7.2f %7.3f" % (range_med, 
                                                                   data_med,
                                                                   r_slope,
                                                                   right_y_int)

# Seventh, we compute root mean square error of the residuals
residuals = [MISS]*n_all # residuals of the fit
fit = [MISS]*n_all       # fitted regression line
valid_count = 0          # total number of non-missing values used
r_sum_sqr_x = 0.0     
r_sum_sqr_e = 0.0        # sum square of residuals
                         # r_slope - slope of linear regression line
                         # r_t - slope error
for i in range(n_all):
    if all_data[i] != MISS:
        valid_count = valid_count + 1
        if valid_count < n_left:
            y_int = left_y_int
        else:
            y_int = right_y_int
        residuals[i] = (y_int + r_slope*(i+1)) - all_data[i]
        fit[i] = y_int + r_slope*(i+1)
        r_sum_sqr_e = r_sum_sqr_e + residuals[i]**2
        r_sum_sqr_x = r_sum_sqr_x + (float(i+1) - data_med)**2

r_se_sqr = r_sum_sqr_e / (valid_count - 2)
r_sb = sqrt(r_se_sqr / r_sum_sqr_x)
r_t = r_slope / r_sb

############## END KENDALLTHIELL()
r_mu = (left_y_int, right_y_int)
r_alpha = r_slope
SSE_red = r_sum_sqr_e

# 3) Now it looks like we compute residuals for all our data
r_residuals = [MISS]*nobs
rssx = 0.0
rsse = [0.0, 0.0]
for k in range(nobs):
    if all_data[k] != MISS:
        if k < (bp+1-left):
            ind = 0
        else:
            ind = 1
        r_residuals[k] = all_data[k] - r_mu[ind] - r_alpha*(k+1)
        rsse[ind] = rsse[ind] + r_residuals[k]**2
        rssx = rssx + (float(k+1) - all_mean)**2

# 4) We now have the squared error and are basically done!
r_sum_sqr_tot = sum(rsse)

############## END KTHTPR1()

# At this point, we have a few things - 
#    r_mu - the y_intercepts of each segment
#    r_alpha - the slope
#    r_sum_sqr_tot - sum sqr total of residuals
#
# We now print out the info and compute critical values, BIC
count = n_all
left_count, right_count = n_left, n_right
sseful = r_sum_sqr_tot # just computed
# F-statistic
f_val = ((sseredmed-sseful)/1.)/(sseful/(count-3))
f_crit = lookup_critical(count-3, "f1")
qslr1, rsq1, rsq2 = bayes(count, sseful, 4)
# amplitude change estimate
y1 = r_mu[0] + r_alpha * range_all[bp+1-left]
y2 = r_mu[1] + r_alpha * range_all[bp+1-right]
est = y1-y2
# k, we have finished this god-awful changepoint
# output string
head = "%7s %6.2f %7.2f %7.2f" % (cmodel, qslr1, rsq1, rsq2)
stats =  " %7.2f %7.2f %7.3f ------- %7.2f %7.2f" % (r_mu[0], r_mu[1],
                                                     r_alpha, f_val, f_crit)
tail = "% 7.2f %5d %4d" % (q_off, n_left, n_right)
print (head+stats+tail)

################################################################################     

## GOOD TO HERE! 


