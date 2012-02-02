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
This module contains representations of various dividing surfaces used in an
RPMD calculation.
"""

import math
import numpy

from rpmd._surface import *
import rpmd.constants as constants

################################################################################

class TransitionState:
    """
    A representation of a transition state dividing surface, defined by a set
    of forming bonds and a set of breaking bonds. The attributes are:
    
    ======================= ====================================================
    Attribute               Description
    ======================= ====================================================
    `geometry`              An array of the coordinates of each atom in each transition state in atomic units
    `formingBonds`          A list of indices for each bond that is forming in the reaction for each transition state
    `breakingBonds`         A list of indices for each bond that is forming in the reaction for each transition state
    ----------------------- ----------------------------------------------------
    `formingBondLengths`    A list of the bond lengths for each bond that is forming in each transition state
    `breakingBondLengths`   A list of the bond lengths for each bond that is breaking in each transition state
    ======================= ====================================================
    
    The `formingBondLengths` and `breakingBondLengths` arrays are automatically
    computed from the given transition state geometries.
    """

    def __init__(self, geometry, formingBonds, breakingBonds):
        
        self.geometry = geometry

        self.formingBonds = numpy.array(formingBonds, numpy.int)
        self.breakingBonds = numpy.array(breakingBonds, numpy.int)
        
        Nts = self.geometry.shape[2]
        Nforming_bonds = formingBonds.shape[1]
        Nbreaking_bonds = breakingBonds.shape[1]
        
        self.formingBondLengths = numpy.empty((Nts, Nforming_bonds))
        self.breakingBondLengths = numpy.empty((Nts, Nbreaking_bonds))
        
        for n in range(Nts):
            for m in range(Nforming_bonds):
                atom1 = self.formingBonds[n,m,0] - 1
                atom2 = self.formingBonds[n,m,1] - 1
                Rx = self.geometry[0,atom1,n] - self.geometry[0,atom2,n]
                Ry = self.geometry[1,atom1,n] - self.geometry[1,atom2,n]
                Rz = self.geometry[2,atom1,n] - self.geometry[2,atom2,n]
                R = math.sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
                self.formingBondLengths[n,m] = R
            
            for m in range(Nbreaking_bonds):
                atom1 = self.breakingBonds[n,m,0] - 1
                atom2 = self.breakingBonds[n,m,1] - 1
                Rx = self.geometry[0,atom1,n] - self.geometry[0,atom2,n]
                Ry = self.geometry[1,atom1,n] - self.geometry[1,atom2,n]
                Rz = self.geometry[2,atom1,n] - self.geometry[2,atom2,n]
                R = math.sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
                self.breakingBondLengths[n,m] = R

    def activate(self):
        """
        Set this object as the active transition state dividing surface in the
        Fortran layer.
        """
        Nts = self.formingBonds.shape[0]
        Nforming_bonds = self.formingBonds.shape[1]
        Nbreaking_bonds = self.breakingBonds.shape[1]

        transition_state.number_of_transition_states = Nts
        transition_state.number_of_forming_bonds = Nforming_bonds
        transition_state.forming_bonds[0:Nts,0:Nforming_bonds,:] = self.formingBonds
        transition_state.forming_bond_lengths[0:Nts,0:Nforming_bonds] = self.formingBondLengths
        transition_state.number_of_breaking_bonds = Nbreaking_bonds
        transition_state.breaking_bonds[0:Nts,0:Nbreaking_bonds,:] = self.breakingBonds
        transition_state.breaking_bond_lengths[0:Nts,0:Nbreaking_bonds] = self.breakingBondLengths
    
    def value(self, position):
        """
        Return the value of the dividing surface function at the given
        `position`.
        """
        self.activate()
        return transition_state.value(position)
    
    def gradient(self, position):
        """
        Return the gradient of the dividing surface function at the given
        `position`.
        """
        self.activate()
        return transition_state.gradient(position)
    
    def hessian(self, position):
        """
        Return the Hessian of the dividing surface function at the given
        `position`.
        """
        self.activate()
        return transition_state.hessian(position)
