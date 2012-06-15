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

import numpy
import unittest

from rpmdrate.interpolate import *

################################################################################

class TestInterpolation(unittest.TestCase):
    """
    Contains unit tests of the interpolation classes.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.xdata = numpy.array([1,2,3,4,5], numpy.float)
        self.ydata = numpy.array([10,23,31,42,49], numpy.float)
        self.N = self.xdata.shape[0]
        self.f1 = LinearInterpolator(self.xdata, self.ydata)
        self.f2 = SemiLogXInterpolator(self.xdata, self.ydata)
        self.f3 = SemiLogYInterpolator(self.xdata, self.ydata)
        self.f4 = LogLogInterpolator(self.xdata, self.ydata)
    
    def test_LinearInterpolator(self):
        """
        Test the LinearInterpolator class.
        """
        x = numpy.arange(1.0, 5.0, 0.01)
        for i in range(x.shape[0]):
            yact = self.f1(x[i])
            for n in range(self.N-1):
                if self.xdata[n] <= x[i] <= self.xdata[n+1]:
                    slope = (self.ydata[n+1] - self.ydata[n]) / (self.xdata[n+1] - self.xdata[n])
                    intercept = self.ydata[n] - slope * self.xdata[n]
                    yexp = slope * x[i] + intercept
                    break
            self.assertAlmostEqual(yexp, yact, 6)

    def test_SemiLogXInterpolator(self):
        """
        Test the SemiLogXInterpolator class.
        """
        x = numpy.arange(1.0, 5.0, 0.01)
        y = numpy.zeros_like(x)
        for i in range(x.shape[0]):
            yact = self.f2(x[i])
            yexp = self.f1(x[i])
            self.assertTrue(0.9 < yexp / yact < 1.1)

    def test_SemiLogYInterpolator(self):
        """
        Test the SemiLogYInterpolator class.
        """
        x = numpy.arange(1.0, 5.0, 0.01)
        y = numpy.zeros_like(x)
        for i in range(x.shape[0]):
            yact = self.f3(x[i])
            yexp = self.f1(x[i])
            self.assertTrue(0.9 < yexp / yact < 1.1)

    def test_LogLogInterpolator(self):
        """
        Test the LogLogInterpolator class.
        """
        x = numpy.arange(1.0, 5.0, 0.01)
        y = numpy.zeros_like(x)
        for i in range(x.shape[0]):
            yact = self.f4(x[i])
            yexp = self.f1(x[i])
            self.assertTrue(0.98 < yexp / yact < 1.02)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
