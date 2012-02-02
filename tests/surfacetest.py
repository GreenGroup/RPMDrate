#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RPMD - Ring polymer molecular dynamics simulations
#
#   Copyright (c) 2012 by Joshua W. Allen (jwallen@mit.edu)
#                         Yury V. Suleimanov (ysuleyma@mit.edu)
#                         William H. Green (whgreen@mit.edu)
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
This script contains unit tests of the :mod:`rpmd.surface` module.
"""

import unittest
import numpy

from rpmd.surface import *

################################################################################

class TestTransitionState(unittest.TestCase):
    """
    Contains unit tests of the :class:`TransitionState` class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        # This is a TS geometry for H + CH4 obtained at the CBS-QB3 level
        self.geometry = numpy.array([
            [[ 0.0000,  0.0000,  0.0000], 
             [ 0.0141, -0.0309,  1.0852], 
             [ 1.0023, -0.0035,  1.5341], 
             [-0.5350,  1.2264,  1.4469], 
             [-0.6769, -0.7373,  1.5351], 
             [-0.8793,  2.0147,  1.6736]],       
        ]).T * 0.52918
        self.formingBonds = numpy.array([
            [(4,6)],
        ], numpy.int)
        self.breakingBonds = numpy.array([
            [(2,4)],
        ], numpy.int)
        self.transitionState = TransitionState(
            geometry = self.geometry,
            formingBonds = self.formingBonds,
            breakingBonds = self.breakingBonds,
        )
        
        # This geometry moves the intermediate H closer to the C and farther
        # from the H (i.e. toward the "reactants" H + CH4)
        self.reactantGeometry = numpy.array([
            [ 0.0000,  0.0000,  0.0000], 
            [ 0.0141, -0.0309,  1.0852], 
            [ 1.0023, -0.0035,  1.5341], 
            [-0.4801,  1.1007,  1.4107], 
            [-0.6769, -0.7373,  1.5351], 
            [-0.8793,  2.0147,  1.6736],       
        ]).T * 0.52918
        # This geometry moves the intermediate H farther from the C and closer
        # to the H (i.e. toward the "products" H2 + CH3)
        self.productGeometry = numpy.array([
            [ 0.0000,  0.0000,  0.0000], 
            [ 0.0141, -0.0309,  1.0852], 
            [ 1.0023, -0.0035,  1.5341], 
            [-0.5899,  1.3521,  1.4831], 
            [-0.6769, -0.7373,  1.5351], 
            [-0.8793,  2.0147,  1.6736],       
        ]).T * 0.52918
    
    def test_geometry(self):
        """
        Test that the TransitionState geometry property was properly set.
        """
        self.assertEqual(self.transitionState.geometry.shape[0], 3)
        self.assertEqual(self.transitionState.geometry.shape[1], 6)
        self.assertEqual(self.transitionState.geometry.shape[2], 1)
        
    def test_value_reactants(self):
        """
        Test the TransitionState.value() method for a geometry in the
        reactant region.
        """
        self.assertTrue(self.transitionState.value(self.reactantGeometry) < 0)

    def test_value_surface(self):
        """
        Test the TransitionState.value() method for a geometry on the
        dividing surface.
        """
        geometry = self.geometry[:,:,0]
        self.assertTrue(abs(self.transitionState.value(geometry)) < 1e-8)

    def test_value_products(self):
        """
        Test the TransitionState.value() method for a geometry in the 
        product region.
        """
        self.assertTrue(self.transitionState.value(self.productGeometry) > 0)

    def test_gradient(self):
        """
        Test the TransitionState.gradient() method by comparing it to the
        numerical gradient.
        """
        dx = 0.0001
        geometry0 = self.geometry[:,:,0]
        gradient = self.transitionState.gradient(geometry0)
        for i in range(3):
            for j in range(6):
                geometry = geometry0.copy()
                geometry[i,j] += dx
                shigh = self.transitionState.value(geometry)
                geometry[i,j] -= 2 * dx
                slow = self.transitionState.value(geometry)
                grad_numerical = (shigh - slow) / (2 * dx)
                grad_analytical = gradient[i,j]
                if grad_analytical == 0:
                    self.assertTrue(abs(grad_numerical) < 1e-6)
                else:
                    self.assertAlmostEqual(grad_analytical, grad_numerical, 6)

    def test_hessian(self):
        """
        Test the TransitionState.hessian() method by comparing it to the
        numerical gradient.
        """
        dx = 0.0002; dy = 0.0001
        geometry0 = self.geometry[:,:,0]
        hessian = self.transitionState.hessian(geometry0)
        for i1 in range(3):
            for j1 in range(6):
                for i2 in range(3):
                    for j2 in range(6):
                        geometry = geometry0.copy()
                        geometry[i1,j1] += dx
                        geometry[i2,j2] += dy
                        shighhigh = self.transitionState.value(geometry)
                        geometry[i1,j1] -= 2 * dx
                        slowhigh = self.transitionState.value(geometry)
                        geometry[i2,j2] -= 2 * dy
                        slowlow = self.transitionState.value(geometry)
                        geometry[i1,j1] += 2 * dx
                        shighlow = self.transitionState.value(geometry)
                        
                        hess_numerical = (shighhigh - shighlow - slowhigh + slowlow) / (4 * dx * dy)
                        hess_analytical = hessian[i1,j1,i2,j2]
                        if hess_analytical == 0:
                            self.assertTrue(abs(hess_numerical) < 1e-6)
                        else:
                            self.assertAlmostEqual(hess_analytical, hess_numerical, 6)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
