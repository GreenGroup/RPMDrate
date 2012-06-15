#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RPMDrate - Bimolecular reaction rates via ring polymer molecular dynamics
#
#   Copyright (c) 2012 by Joshua W. Allen (jwallen@mit.edu)
#                         William H. Green (whgreen@mit.edu)
#                         Yury V. Suleimanov (ysuleyma@mit.edu, ysuleyma@princeton.edu)
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
This script contains unit tests of the :mod:`rpmd.constants` module.
"""

import unittest
import math

import rpmdrate.constants as constants

################################################################################

class TestConstants(unittest.TestCase):
    """
    Contains unit tests that ensure that the physical constants are visible in
    both pure Python and compiled Cython modes, and in both cases are set to
    the appropriate values.
    """
    
    def test_avogadroConstant(self):
        """
        Test the value of the Avogadro constant.
        """
        Na = 6.02214179e23
        self.assertAlmostEqual(constants.Na / Na, 1.0, 6, '{0} != {1}'.format(constants.Na, Na))

    def test_boltzmannConstant(self):
        """
        Test the value of the Boltzmann constant.
        """
        kB = 1.3806504e-23
        self.assertAlmostEqual(constants.kB / kB, 1.0, 6, '{0} != {1}'.format(constants.kB, kB))
    
    def test_elementaryCharge(self):
        """
        Test the value of the elementary charge constant.
        """
        e = 1.602176565e-19
        self.assertAlmostEqual(constants.e / e, 1.0, 6, '{0} != {1}'.format(constants.e, e))
    
    def test_gasLawConstant(self):
        """
        Test the value of the gas law constant.
        """
        R = 8.314472
        self.assertAlmostEqual(constants.R / R, 1.0, 6, '{0} != {1}'.format(constants.R, R))
    
    def test_planckConstant(self):
        """
        Test the value of the Planck constant.
        """
        h = 6.62606896e-34
        self.assertAlmostEqual(constants.h / h, 1.0, 6, '{0} != {1}'.format(constants.h, h))
    
    def test_reducedPlanckConstant(self):
        """
        Test the value of the reduced Planck constant.
        """
        hbar = 1.054571726e-34
        self.assertAlmostEqual(constants.hbar / hbar, 1.0, 6, '{0} != {1}'.format(constants.hbar, hbar))
    
    def test_pi(self):
        """
        Test the value of pi.
        """
        self.assertAlmostEqual(constants.pi / math.pi, 1.0, 6, '{0} != {1}'.format(constants.pi, math.pi))

    def test_speedOfLight(self):
        """
        Test the value of the speed of light in a vacuum.
        """
        c = 299792458
        self.assertEqual(constants.c, c)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
