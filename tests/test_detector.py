# -*- coding: utf-8 -*-
"""
Unit tests for the radio_simus.detector module

Usage: python3.7 tests/test_detector.py


"""

import unittest
import sys
import numpy as np

import astropy.units as u

from os.path import split, join, realpath
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))
from radio_simus.detector import Detector


class DetectorTest(unittest.TestCase):
    """Unit tests for the detector module

    TODO: add tests for missing modules

    """

    def assertQuantity(self, x, y):
        """Check that two astropy.Quantities are consistent"""
        self.assertEqual((x / x.unit).value, (y / x.unit).value)

    def test_type(self):
        print("Start _type test")

        type_array = np.array([ "LDPA", "butterfly", "Horizon" ])

        det = Detector()
        det.type = type_array[:].tolist()
        types = det.type

        self.assertEqual(type_array[1], types[1]) #shall be equal
        self.assertEqual(len(type_array), len(types))

        print("Finish _type test")


    def test_position(self):
        print("Start _position test")

        pos_array = np.array([ [100, 100, 3000], [200, 300, 2500], [10, 500, -100] ])

        det = Detector()
        det.position = pos_array[:].tolist() # XXX should be an astropy.Quantity
        positions = det.position

        self.assertQuantity(pos_array[1,1] * u.m, positions[1,1]) #shall be equal
        self.assertEqual(len(pos_array), len(positions))

        print("Finish _position test")


    def test_ID(self):
        print("Start _ID test")

        ID_array = np.array([1,43,600,"ant1", 99])

        det = Detector()
        det.ID = ID_array[:].tolist()
        ID = det.ID

        self.assertEqual(ID_array[1], ID[1]) #shall be equal
        self.assertEqual(len(ID_array), len(ID))

        print("Finish _ID test")


    def test_location(self):
        print("Start _location test")
        det = Detector()

        det.location = "Fancy place"

        self.assertEqual("Fancy place", det.location) #shall be equal

        print("Finish _location test")


    def test_origin(self):
        print("Start _origin test")
        det = Detector()

        origin_default = ( 0 * u.m, 0 * u.m, 0 * u.m)
        det.origin = origin_default

        self.assertQuantity(origin_default[0], det.origin[0]) #shall be equal

        print("Finish _origin test")


if __name__ == "__main__":
    unittest.main()
