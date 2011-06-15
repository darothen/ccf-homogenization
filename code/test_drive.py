import ushcn_io
import random
from operator import itemgetter

from ushcn_data import Station, Series
from util import compute_arc_dist, compute_monthly_anomalies
from util import compute_first_diff, compute_corr, compute_std

data_src = "raw"
variable = "max"

nstns = 100
hidist = 200.
distinc = 200.
numsrt = 20
begyr = 1920
endyr = 1980

all_series, all_stations = ushcn_io.get_ushcn_data(data_src, variable)

station_ids = sorted(random.sample(all_stations.keys(), nstns))
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

    sorted_neighbors = sorted(neighbor_dict.iteritems(), key=itemgetter(1))
    
    close_neighbors = []
    iter = 1
    mindist = hidist
    print "...searching for neighbors within %4.1fkm" % mindist
    while (len(close_neighbors) < numsrt):
        print "......iteration %d, %d neighbors found" % (iter, 
                                                          len(close_neighbors))
        for (coop_id, dist) in sorted_neighbors:
            if dist < mindist:
                sorted_neighbors.remove((coop_id, dist))
                close_neighbors.append((coop_id, dist))
        mindist = mindist+distinc
        print ".........distance now %4.1fkm" % mindist    
        iter = iter+1
    print "......ultimately found %d neighbors" % len(close_neighbors)
    all_neighbors[coop_id1] = close_neighbors

##########################################################################
years = range(begyr, endyr)
for s in series.itervalues():
    trunc_series = s.trunc_series(begyr, endyr)
    flat_trunc = []
    for monthly in trunc_series:
        flat_trunc.extend(monthly[:-1])
    anomalies = compute_monthly_anomalies(flat_trunc, -9999)
    s.set_series(anomalies, years)
    
print "Determining correlated neighbors"

##### WARNING WARNING WARNING #####
## It still might be possible that the truncated series are **too short**. 
## Will need to go into ushchn_data and be sure that the data is padded to fill
## out begyr to endyr.

all_corrs = dict()
for coop_id1 in station_ids:
    
    print "...%s" % coop_id1
    neighbors = [n[0] for n in all_neighbors[coop_id1]]
    for coop_id2 in neighbors:
        print "......%s" % coop_id2
        
        cand = series[coop_id1]
        cand_data = cand.monthly_series
        neighb = series[coop_id2]
        neighb_data = neighb.monthly_series
        
        MISS = cand.MISSING_VAL
        
        print ".........Aligning cand/neighb series"
        cand_align, neighb_align = [], []
        for (cand_val, neighb_val) in zip(cand_data, neighb_data):
            if (cand_val != MISS and neighb_val != MISS):
                cand_align.append(cand_val)
                neighb_align.append(neighb_val)
        print "            %d %d" % (len(cand_align), len(neighb_align))
                
        print ".........Computing first differences"
        cand_dif = compute_first_diff(cand_align, MISS)
        neighb_dif = compute_first_diff(neighb_align, MISS)
        
        print ".........Computing correlation coefficient"
        r = compute_corr(cand_dif, neighb_dif, MISS)
        cand_std = compute_std(cand_dif, MISS)
        neighb_std = compute_std(neighb_dif, MISS)
        print "            %1.3f %3.3f %3.3f" % (r, cand_std, neighb_std)

network = dict(stations=stations, series=series, neighbors=all_neighbors)

