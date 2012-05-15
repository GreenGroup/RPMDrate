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
This module contains functions for unit conversions of various types of
physical quantities.
"""

import numpy

import rpmdrate.constants as constants

################################################################################

TEMPERATURE_UNITS = {
    'K': 1.0,
    'degC': 1.0,
    'degF': 5.0/9.0,
    'degR': 5.0/9.0,
}

def convertTemperature(quantity, units):
    """
    Convert a given `quantity` with units of temperature to the given `units` 
    of temperature. A :class:`ValueError` is raised if this conversion is not
    successful.
    """
    if isinstance(quantity, tuple):
        value0, units0 = quantity
    else:
        raise ValueError('Invalid value "{0}" for quantity; must be a tuple (value,units) with units of temperature.'.format(quantity))
    
    if isinstance(value0, (list, tuple)):
        value0 = numpy.array(value0)
    
    if units0 not in TEMPERATURE_UNITS:
        raise ValueError('Invalid input units "{0}" for temperature.'.format(units0))
    if units not in TEMPERATURE_UNITS:
        raise ValueError('Invalid output units "{0}" for temperature.'.format(units))
    
    if units0 == 'degF':        value0 += 459.67
    elif units0 == 'degC':      value0 += 273.15
    
    value = value0 * TEMPERATURE_UNITS[units0] / TEMPERATURE_UNITS[units]
    
    if units == 'degF':         value -= 459.67
    elif units == 'degC':       value -= 273.15
     
    return value

################################################################################

TIME_UNITS = {
    'fs': 1.0e-15,
    'ps': 1.0e-12,
    'ns': 1.0e-9,
    'us': 1.0e-6,
    'ms': 1.0e-3,
    's': 1.0,
    'min': 60.0,
    'hr': 3600.0,
    'day': 86400.0,
}

def convertTime(quantity, units):
    """
    Convert a given `quantity` with units of time to the given `units` 
    of time. A :class:`ValueError` is raised if this conversion is not
    successful.
    """
    if isinstance(quantity, tuple):
        value0, units0 = quantity
    else:
        raise ValueError('Invalid value "{0}" for quantity; must be a tuple (value,units) with units of time.'.format(quantity))
    
    if isinstance(value0, (list, tuple)):
        value0 = numpy.array(value0)
    
    if units0 not in TIME_UNITS:
        raise ValueError('Invalid input units "{0}" for time.'.format(units0))
    if units not in TIME_UNITS:
        raise ValueError('Invalid output units "{0}" for time.'.format(units))
    
    return value0 * TIME_UNITS[units0] / TIME_UNITS[units]

################################################################################

LENGTH_UNITS = {
    'fm': 1.0e-15,
    'pm': 1.0e-12,
    'nm': 1.0e-9,
    'um': 1.0e-6,
    'mm': 1.0e-3,
    'cm': 1.0e-2,
    'dm': 1.0e-1,
    'm': 1.0,
    'km': 1.0e3,
    'angstrom': 1.0e-10,
    'bohr': 5.2917721092e-11,
}

def convertLength(quantity, units):
    """
    Convert a given `quantity` with units of length to the given `units` 
    of length. A :class:`ValueError` is raised if this conversion is not
    successful.
    """
    if isinstance(quantity, tuple):
        value0, units0 = quantity
    else:
        raise ValueError('Invalid value "{0}" for quantity; must be a tuple (value,units) with units of length.'.format(quantity))
    
    if isinstance(value0, (list, tuple)):
        value0 = numpy.array(value0)
    
    if units0 not in LENGTH_UNITS:
        raise ValueError('Invalid input units "{0}" for length.'.format(units0))
    if units not in LENGTH_UNITS:
        raise ValueError('Invalid output units "{0}" for length.'.format(units))
    
    return value0 * LENGTH_UNITS[units0] / LENGTH_UNITS[units]

################################################################################

MASS_UNITS = {
    'mg': 1.0e-6,
    'g': 1.0e-3,
    'kg': 1.0,
}
MOLAR_MASS_UNITS = {
    'mg/mol': 1.0e-6 / constants.Na,
    'g/mol': 1.0e-3 / constants.Na,
    'kg/mol': 1.0 / constants.Na,
    'mg/kmol': 1.0e-9 / constants.Na,
    'g/kmol': 1.0e-6 / constants.Na,
    'kg/kmol': 1.0e-3 / constants.Na,
    'amu': 1.0e-3 / constants.Na,
}

def convertMass(quantity, units):
    """
    Convert a given `quantity` with units of mass to the given `units` 
    of mass. A :class:`ValueError` is raised if this conversion is not
    successful.
    """
    if isinstance(quantity, tuple):
        value0, units0 = quantity
    else:
        raise ValueError('Invalid value "{0}" for quantity; must be a tuple (value,units) with units of mass.'.format(quantity))
    
    if isinstance(value0, (list, tuple)):
        value0 = numpy.array(value0)
    
    if units0 not in MASS_UNITS and units0 not in MOLAR_MASS_UNITS:
        raise ValueError('Invalid input units "{0}" for mass.'.format(units0))
    if units not in MASS_UNITS and units not in MOLAR_MASS_UNITS:
        raise ValueError('Invalid output units "{0}" for mass.'.format(units))
    
    if units0 in MASS_UNITS and units in MASS_UNITS:
        return value0 * MASS_UNITS[units0] / MASS_UNITS[units]
    elif units0 in MASS_UNITS and units in MOLAR_MASS_UNITS:
        return value0 * MASS_UNITS[units0] / MOLAR_MASS_UNITS[units]
    elif units0 in MOLAR_MASS_UNITS and units in MASS_UNITS:
        return value0 * MOLAR_MASS_UNITS[units0] / MASS_UNITS[units]
    elif units0 in MOLAR_MASS_UNITS and units in MOLAR_MASS_UNITS:
        return value0 * MOLAR_MASS_UNITS[units0] / MOLAR_MASS_UNITS[units]
    
    return value0 * LENGTH_UNITS[units0] / LENGTH_UNITS[units]
