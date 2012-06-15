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
This module contains representations of various thermostats available in RPMD.
"""

import os.path
import numpy

import rpmdrate.constants as constants
import rpmdrate.quantity as quantity

################################################################################

class AndersenThermostat:
    """
    A representation of an Andersen thermostat, a simple method of sampling
    the NVT ensemble which periodically replaces the momenta with a fresh
    sampling from a Gaussian distribution at the temperature of interest, as
    if resulting from a collision with a heat bath.
    """
    
    def __init__(self, samplingTime=None):
        if samplingTime is not None:
            self.samplingTime = float(quantity.convertTime(samplingTime, "ps") / 2.418884326505e-5)
        else:
            self.samplingTime = 0.0

    def activate(self, module, Natoms, Nbeads):
        """
        Set the thermostat as active in the Fortran layer of the given
        `module`.
        """
        module.thermostat = 1
        module.andersen_sampling_time = self.samplingTime

################################################################################

class GLEThermostat(object):
    """
    A representation of a colored-noise, generalized Langevin equation
    thermostat. The GLE thermostat offers significantly faster convergence
    when compared to the simpler Andersen thermostat, but requires more effort
    to set up. In particular, the thermostat requires two matrices be specified
    as input; these can be generated at 
    
        http://gle4md.berlios.de/compose.php?page=matrix
    
    and supplied either as a file on disk or as numpy arrays. The A and C
    matrices are assumed to have units of s^-1 and K, respectively.
    """
    
    def __init__(self, A=None, C=None):
        self.A = A
        self.C = C

    @property
    def A(self):
        return self._A
    @A.setter
    def A(self, value):
        if isinstance(value, (list,tuple)) and len(value) == 2 and os.path.exists(value[0]):
            # value is the path of a file on disk to load from and the corresponding units
            path, units = value
            _A = []
            f = open(path, 'r')
            for line in f:
                # Remove comment
                if '#' in line: line = line[0:line.index('#')].strip()
                tokens = line.split()
                if len(tokens) > 0:
                    _A.append([float(t) for t in tokens])
            f.close()
            self._A = numpy.array(quantity.convertFrequency((_A,units),'s^-1'))
        elif isinstance(value, numpy.ndarray):
            self._A = value
        elif value is None:
            self._A = None
        else:
            raise ValueError('Unexpected value {0!r} for A attribute.'.format(value))
        
    @property
    def C(self):
        return self._C
    @C.setter
    def C(self, value):
        if isinstance(value, (list,tuple)) and len(value) == 2 and os.path.exists(value[0]):
            # value is the path of a file on disk to load from and the corresponding units
            path, units = value
            _C = []
            f = open(path, 'r')
            for line in f:
                # Remove comment
                if '#' in line: line = line[0:line.index('#')].strip()
                tokens = line.split()
                if len(tokens) > 0:
                    _C.append([float(t) for t in tokens])
            f.close()
            self._C = numpy.array(quantity.convertTemperature((_C,units),'K'))
        elif isinstance(value, numpy.ndarray):
            self._C = value
        elif value is None:
            self._C = None
        else:
            raise ValueError('Unexpected value {0!r} for C attribute.'.format(value))

    def activate(self, module, Natoms, Nbeads):
        """
        Set the thermostat as active in the Fortran layer of the given
        `module`.
        """
        Ns = self._A.shape[0] - 1
        module.thermostat = 2
        module.gle_ns = Ns
        module.gle_a[0:Ns+1,0:Ns+1] = self._A * 2.418884326505e-17  # s^-1 to atomic units of inverse time
        if self._C is None:
            module.gle_c[0:Ns+1,0:Ns+1] = numpy.zeros((Ns+1,Ns+1))
            for s in range(Ns+1):
                module.gle_c[s,s] = Nbeads / module.beta
        else:
            module.gle_c[0:Ns+1,0:Ns+1] = self._C * constants.kB / 4.35974417e-18  # K to atomic units of energy
