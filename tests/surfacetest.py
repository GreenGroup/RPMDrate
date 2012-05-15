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

from rpmdrate.surface import *

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
        self.geometry = numpy.array(
            [[ 0.0000,  0.0000,  0.0000], 
             [ 0.0141, -0.0309,  1.0852], 
             [ 1.0023, -0.0035,  1.5341], 
             [-0.5350,  1.2264,  1.4469], 
             [-0.6769, -0.7373,  1.5351], 
             [-0.8793,  2.0147,  1.6736]], order='F')
        self.formingBonds = numpy.array([(4,6)], numpy.int)
        self.breakingBonds = numpy.array([(2,4)], numpy.int)
        self.transitionState = TransitionState(
            geometry = (self.geometry,"angstrom"),
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
        ], order='F').T / 0.52918
        # This geometry moves the intermediate H farther from the C and closer
        # to the H (i.e. toward the "products" H2 + CH3)
        self.productGeometry = numpy.array([
            [ 0.0000,  0.0000,  0.0000], 
            [ 0.0141, -0.0309,  1.0852], 
            [ 1.0023, -0.0035,  1.5341], 
            [-0.5899,  1.3521,  1.4831], 
            [-0.6769, -0.7373,  1.5351], 
            [-0.8793,  2.0147,  1.6736],       
        ], order='F').T / 0.52918
    
    def test_geometry(self):
        """
        Test that the TransitionState geometry property was properly set.
        """
        self.assertEqual(self.transitionState.geometry.shape[0], 3)
        self.assertEqual(self.transitionState.geometry.shape[1], 6)
        
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
        geometry = self.geometry.T / 0.52918
        self.assertTrue(abs(self.transitionState.value(geometry)) < 1e-4)

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
        geometry0 = self.geometry.T / 0.52918
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
        geometry0 = self.geometry.T / 0.52918
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

class TestReactants(unittest.TestCase):
    """
    Contains unit tests of the :class:`Reactants` class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.reactants = Reactants(
            atoms = ['H', 'H', 'H'],
            reactant1Atoms = [1],
            reactant2Atoms = [2,3],
            Rinf = (5,"angstrom"),
        )
        self.position = numpy.array([
            [ 0.02,  0.04,  -0.929764], 
            [ -0.03,  0.01,  0.07], 
            [ 0.05,  0.0000,  0.929764], 
        ]).T / 0.52918
    
    def test_value(self):
        """
        Test the Reactants.value() method for a geometry.
        """
        self.assertAlmostEqual(self.reactants.value(self.position), 6.74608, 4)

    def test_gradient(self):
        """
        Test the Reactants.gradient() method by comparing it to the
        numerical gradient.
        """
        dx = 0.0001
        
        ds0 = self.reactants.gradient(self.position)
        for i in range(3):
            for j in range(3):
                position = self.position.copy()
                position[i,j] += dx
                s0_high = self.reactants.value(position)
                position[i,j] -= 2 * dx
                s0_low = self.reactants.value(position)

                ds0_exp = (s0_high - s0_low) / (2 * dx)
                ds0_act = ds0[i,j]
                self.assertAlmostEqual(ds0_exp / ds0_act, 1.0, 6, '{0} != {1}'.format(ds0_exp, ds0_act))

    def test_hessian(self):
        """
        Test the Reactants.hessian() method by comparing it to the
        numerical Hessian.
        """
        dx = 0.002; dy = 0.001
        
        d2s0 = self.reactants.hessian(self.position)
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    for d in range(3):
                        position = self.position.copy()
                        position[a,b] += dx; position[c,d] += dy
                        s0_hh = self.reactants.value(position)
                        position[c,d] -= 2 * dy
                        s0_hl= self.reactants.value(position)
                        position[a,b] -= 2 * dx
                        s0_ll = self.reactants.value(position)
                        position[c,d] += 2 * dy
                        s0_lh = self.reactants.value(position)
                        
                        d2s0_exp = (s0_hh - s0_hl - s0_lh + s0_ll) / (4 * dx * dy)
                        d2s0_act = d2s0[a,b,c,d]
                        self.assertAlmostEqual(d2s0_exp / d2s0_act, 1.0, 4, '{0} != {1}'.format(d2s0_exp, d2s0_act))

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
