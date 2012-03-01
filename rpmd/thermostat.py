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
This module contains representations of various thermostats available in RPMD.
"""

import rpmd.quantity as quantity

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

    def activate(self, module):
        """
        Set the thermostat as active in the Fortran layer of the given
        `module`.
        """
        module.thermostat = 1
        module.andersen_sampling_time = self.samplingTime
