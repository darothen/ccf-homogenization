#!/usr/bin/env python
#
# Daniel Rothenberg, 2011-06-14

"""Test cases for straight-forward mathematical computations 
involved in station network analysis.

"""
__docformat__ = "restructuredtext"

# http://docs.python.org/library/unittest.html
import unittest
# http://docs.python.org/library/math.html
from math import pi
# ccf-homogenization imports
from ushcn_data import Station
from parameters import RADIUS_EARTH
from util import compute_arc_dist, compute_mean, anomaly
from util import compute_corr, compute_std

class TimeseriesMathChecks(unittest.TestCase):
    
    #: These are test monthly datasets for computing anomalies
    dataset1 = [20, 30, -9999, 50, 60]
    dataset2 = [20, 30, 40, 50, 60, 70]
    dataset3 = [-9999, -9999, -9999, -9999]
    
    def testDataset1(self):
        """Dataset with one missing value"""        
        mean = compute_mean(self.dataset1, -9999)
        self.assertEquals(mean, 40.0)
        
        std = compute_std(self.dataset1, -9999)
        self.assertAlmostEquals(std, 18.2574, delta=1e-3)
        
    def testDataset2(self):
        """Dataset with no missing values"""        
        mean = compute_mean(self.dataset2, -9999)
        self.assertEquals(mean, 45.0)
        
        std = compute_std(self.dataset2, -9999)
        self.assertAlmostEquals(std, 18.7082, delta=1e-3)
        
    def testDataset3(self):
        """Dataset with only missing values"""        
        mean = compute_mean(self.dataset3, -9999)
        self.assertEquals(mean, -9999)
        
        std = compute_mean(self.dataset3, -9999)
        self.assertAlmostEquals(std, -9999, delta=1e-3)
        
    def testCorrelation1(self):
        """A random dataset found at 
        http://www.stat.wmich.edu/s216/book/node122.html"""
        x = [1, 3, 4, 4]
        y = [2, 5, 5, 8]
        r = compute_corr(x, y, valid=True)
        self.assertAlmostEquals(r, 0.866, delta=1e-3)
        
    def testCorrelation2(self):
        """Another random dataset from 
        http://www.socialresearchmethods.net/kb/statcorr.php"""
        x = [68, 71, 62, 75, 58, 60, 67, 68, 71, 69, 68, 67, 63, 62, 60, 63, 65,
             67, 63, 61]
        y = [4.1, 4.6, 3.8, 4.4, 3.2, 3.1, 3.8, 4.1, 4.3, 3.7, 3.5, 3.2, 3.7, 3.3,
             3.4, 4.0, 4.1, 3.8, 3.4, 3.6]
        r = compute_corr(x, y, valid=True)
        self.assertAlmostEquals(r, 0.73, delta=1e-3)
        
    def testCorrelation3(self):
        """Too little data should yield correlation of None"""
        x = [1, ]
        y = [2, ]
        r = compute_corr(x, y, -9999)
        self.assertIsNone(r)
        
    def testCorrelation4(self):
        """Correlation between same dataset should be 1.0"""
        x = [1-1e9, 0, 1+1e9]
        r = compute_corr(x, x, valid=True)
        self.assertAlmostEquals(r, 1.0, delta=1e-3)
    
    def testCorrelation5(self):
        """Correlation between datasets with standard deviation
        of 0 should raise a ZeroDivisonError"""
        x = [1-1e9, 1-1e9, 1-1e9]
        y = [1+1e9, 1+1e9, 1+1e9]
        self.assertRaises(ZeroDivisionError, compute_corr, x, y, valid=True)
        
    def testCorrelation6(self):
        x1, x2 = [1,2,-9999,4,5], [1,4,5]
        y1, y2 = [1,-9999,3,5,4], [1,5,4]
        r1 = compute_corr(x1,y1,-9999)
        r2 = compute_corr(x2,y2,-9999)
        self.assertEquals(r1,r2)
        
    def testStd1(self):
        """Shouldn't be able to compute std if less than 2 values"""
        miss = -9999
        data = [3.3]
        std = compute_std(data, miss, valid=True)
        self.assertEquals(std, miss)
        
    def testStd2(self):
        """Should return the missing value if empty dataset given."""
        miss = -9999
        data = []
        std = compute_std(data, miss)
        self.assertEquals(std, miss)
        
    def anomalies(self):
        """Compute anomalies on a valid datum and invalid datum."""
        x1 = 43
        x2 = -9999
        mean = 46
        self.assertEquals(-3, anomaly(x1, mean))
        self.assertEquals(-9999, anomaly(x2, mean))
    

class SphereMathChecks(unittest.TestCase):
    
    #: These are stations to test distances between.
    lat_a, lon_a = (31.0581, -87.0547)
    STATION_A = Station(lat=lat_a, lon=lon_a) # BREWTON, AL
    lat_b, lon_b = (42.45, -76.45)
    STATION_B = Station(lat=lat_b, lon=lon_b) # ITHACA, NY
    lat_c, lon_c = 90.0, 0.1
    STATION_C = Station(lat=lat_c, lon=lon_c) # FAKE NORTH POLE STATION
    lat_d, lon_d = 0., 0.
    STATION_D = Station(lat=lat_d, lon=lon_d) # FAKE EQUATOR STATION
    
    #: The delta tolerance between calculated and expected distances.
    #: This is very conservative because the values found online may
    #: be generated using a different mathematical formula which
    #: may not be based on a spherical Earth or may be contaminated
    #: with some round-off error.
    DELTA = 5.
    
    def testSameStation(self):
        """Distance between the same station should round off to 0."""
        dist = compute_arc_dist(self.STATION_A, self.STATION_A)
        self.assertAlmostEquals(dist, 0., delta=self.DELTA)
        
    def testTwoStations(self):
        """According to jan.ucc.nau.edu/~cvm/latlongdist, the distance
        between Ithaca and Brewton is 1576.7066 km
        """
        dist1 = compute_arc_dist(self.STATION_A, self.STATION_B)
        dist2 = compute_arc_dist(self.STATION_B, self.STATION_A)
        self.assertAlmostEquals(dist1, 1576.7066, delta=self.DELTA)
        self.assertAlmostEquals(dist2, 1576.7066, delta=self.DELTA)
        
    def testQuarterArc(self):
        dist1 = compute_arc_dist(self.STATION_C, self.STATION_D)
        dist2 = compute_arc_dist(self.STATION_D, self.STATION_C)
        quarter_arc = (2.*pi*RADIUS_EARTH)/4.0
        self.assertAlmostEquals(dist1, quarter_arc, delta=self.DELTA)
        self.assertAlmostEquals(dist2, quarter_arc, delta=self.DELTA)