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
from util import compute_arc_dist, compute_monthly_anomalies

class TimeseriesChecks(unittest.TestCase):
    
    #: These are test monthly datasets
    dataset1 = [20, 30, -9999, 50, 60]
    dataset2 = [20, 30, 40, 50, 60, 70]
    dataset3 = [-9999, -9999, -9999, -9999]
    
    def testDataset1(self):
        """Dataset with one missing value"""
        anomalies = compute_monthly_anomalies(self.dataset1, -9999.)
        self.assertEquals(anomalies, [-20, -10, -9999, 10, 20])
        
    def testDataset2(self):
        """Dataset with no missing values"""
        anomalies = compute_monthly_anomalies(self.dataset2, -9999.)
        self.assertEquals(anomalies, [-25, -15, -5, 5, 15, 25])
        
    def testDataset3(self):
        """Dataset with only missing values"""
        anomalies = compute_monthly_anomalies(self.dataset3, -9999.)
        self.assertEquals(anomalies, [-9999, -9999, -9999, -9999])

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