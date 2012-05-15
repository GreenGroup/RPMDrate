#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RPMDrate - Ring polymer molecular dynamics simulations
#
#   Copyright (c) 2012 by Joshua W. Allen (jwallen@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a 
#   copy of this software and associated documentation files (the "Software"), 
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
#   and/or sell copies of the Software, and to permit persons to whom the 
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
#   DEALINGS IN THE SOFTWARE. 
#
################################################################################

"""
This script contains unit tests of the :mod:`rpmdrate.quantity` module.
"""

import unittest
import math
import quantities as pq

from rpmdrate.quantity import *

################################################################################

class TestConvertLength(unittest.TestCase):
    """
    Contains unit tests of the :func:`convertLength` method.
    """

    def test_m_to_cm(self):
        """
        Test the conversion of a length from m to cm.
        """
        self.assertAlmostEqual(convertLength((1.0,"m"),"cm"), 1.0e2, 6)

    def test_cm_to_m(self):
        """
        Test the conversion of a length from cm to m.
        """
        self.assertAlmostEqual(convertLength((1.0,"cm"),"m"), 1.0e-2, 6)

    def test_m_to_angstrom(self):
        """
        Test the conversion of a length from m to cm.
        """
        self.assertAlmostEqual(convertLength((1.0e-10,"m"),"angstrom"), 1.0, 6)

    def test_angstrom_to_m(self):
        """
        Test the conversion of a length from cm to m.
        """
        self.assertAlmostEqual(convertLength((1.0e10,"angstrom"),"m"), 1.0, 6)

################################################################################

class TestConvertMass(unittest.TestCase):
    """
    Contains unit tests of the :func:`convertMass` method.
    """

    def test_g_to_gpermol(self):
        """
        Test the conversion of a mass from g to g/mol.
        """
        self.assertAlmostEqual(convertMass((1.0/constants.Na,"g"),"g/mol"), 1.0, 3)

    def test_gpermol_to_g(self):
        """
        Test the conversion of a mass from g/mol to g.
        """
        self.assertAlmostEqual(convertMass((1.0*constants.Na,"g/mol"),"g"), 1.0, 3)

    def test_g_to_kgperkmol(self):
        """
        Test the conversion of a mass from g to g/mol.
        """
        self.assertAlmostEqual(convertMass((1.0/constants.Na,"g"),"kg/kmol"), 1.0, 3)

    def test_kgperkmol_to_g(self):
        """
        Test the conversion of a mass from g/mol to g.
        """
        self.assertAlmostEqual(convertMass((1.0*constants.Na,"kg/kmol"),"g"), 1.0, 3)

    def test_g_to_amu(self):
        """
        Test the conversion of a mass from g to amu.
        """
        self.assertAlmostEqual(convertMass((1.0/constants.Na,"g"),"amu"), 1.0, 3)

    def test_amu_to_g(self):
        """
        Test the conversion of a mass from amu to g.
        """
        self.assertAlmostEqual(convertMass((1.0*constants.Na,"amu"),"g"), 1.0, 3)

    def test_g_to_kg(self):
        """
        Test the conversion of a mass from g to kg.
        """
        self.assertAlmostEqual(convertMass((1000.,"g"),"kg"), 1.0, 3)

    def test_kg_to_g(self):
        """
        Test the conversion of a mass from kg to g.
        """
        self.assertAlmostEqual(convertMass((0.001,"kg"),"g"), 1.0, 3)

################################################################################

class TestConvertTemperature(unittest.TestCase):
    """
    Contains unit tests of the :func:`convertTemperature` method.
    """
    
    def test_degF_to_degC(self):
        """
        Test the conversion of a temperature from deg F to deg C.
        """
        self.assertAlmostEqual(convertTemperature((68,"degF"),"degC"), 20, 6)

    def test_degF_to_degR(self):
        """
        Test the conversion of a temperature from deg F to deg R.
        """
        self.assertAlmostEqual(convertTemperature((68,"degF"),"degR"), 527.67, 6)

    def test_degF_to_K(self):
        """
        Test the conversion of a temperature from deg F to K.
        """
        self.assertAlmostEqual(convertTemperature((68,"degF"),"K"), 293.15, 6)

    def test_degC_to_degF(self):
        """
        Test the conversion of a temperature from deg C to deg F.
        """
        self.assertAlmostEqual(convertTemperature((20,"degC"),"degF"), 68, 6)

    def test_degC_to_degR(self):
        """
        Test the conversion of a temperature from deg C to deg R.
        """
        self.assertAlmostEqual(convertTemperature((20,"degC"),"degR"), 527.67, 6)

    def test_degC_to_K(self):
        """
        Test the conversion of a temperature from deg C to K.
        """
        self.assertAlmostEqual(convertTemperature((20,"degC"),"K"), 293.15, 6)

    def test_degR_to_degF(self):
        """
        Test the conversion of a temperature from deg R to deg F.
        """
        self.assertAlmostEqual(convertTemperature((527.67,"degR"),"degF"), 68, 6)

    def test_degR_to_degC(self):
        """
        Test the conversion of a temperature from deg R to deg C.
        """
        self.assertAlmostEqual(convertTemperature((527.67,"degR"),"degC"), 20, 6)

    def test_degR_to_K(self):
        """
        Test the conversion of a temperature from deg R to K.
        """
        self.assertAlmostEqual(convertTemperature((527.67,"degR"),"K"), 293.15, 6)

    def test_K_to_degF(self):
        """
        Test the conversion of a temperature from K to deg F.
        """
        self.assertAlmostEqual(convertTemperature((293.15,"K"),"degF"), 68, 6)

    def test_K_to_degC(self):
        """
        Test the conversion of a temperature from K to deg C.
        """
        self.assertAlmostEqual(convertTemperature((293.15,"K"),"degC"), 20, 6)

    def test_K_to_degR(self):
        """
        Test the conversion of a temperature from K to deg K.
        """
        self.assertAlmostEqual(convertTemperature((293.15,"K"),"degR"), 527.67, 6)

################################################################################

class TestConvertTime(unittest.TestCase):
    """
    Contains unit tests of the :func:`convertTime` method.
    """

    def test_s_to_ms(self):
        """
        Test the conversion of a time from s to ms.
        """
        self.assertAlmostEqual(convertTime((0.001,"s"),"ms"), 1.0, 6)

    def test_ms_to_s(self):
        """
        Test the conversion of a time from ms to s.
        """
        self.assertAlmostEqual(convertTime((1000.,"ms"),"s"), 1.0, 6)

    def test_s_to_hr(self):
        """
        Test the conversion of a time from s to hr.
        """
        self.assertAlmostEqual(convertTime((3600.,"s"),"hr"), 1.0, 6)

    def test_hr_to_s(self):
        """
        Test the conversion of a time from hr to s.
        """
        self.assertAlmostEqual(convertTime((1.0,"hr"),"s"), 3600., 6)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
