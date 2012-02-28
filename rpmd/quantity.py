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
physical quantities. The unit conversions make use of the :mod:`quantities`
package.
"""

import quantities as pq

# Define custom units
pq.UnitQuantity('bohr', pq.constants.atomic_unit_of_length, 'bohr')

################################################################################

TEMPERATURE_DIMENSIONS = [pq.K.simplified.dimensionality]

def convertTemperature(quantity, units):
    """
    Convert a given `quantity` with units of temperature to the given `units` 
    of temperature. A :class:`ValueError` is raised if this conversion is not
    successful.
    """
    if isinstance(quantity, tuple):
        quantity = pq.Quantity(quantity[0], quantity[1])
    elif not isinstance(quantity, pq.Quantity):
        raise ValueError('Invalid value "{0}" for quantity; must be a tuple or Quantity object with units of temperature.'.format(quantity))
    
    if isinstance(units, str):
        units = pq.Quantity(1.0, units)
    
    inputDimensionality = quantity.units.dimensionality
    outputDimensionality = units.dimensionality

    if inputDimensionality.simplified not in TEMPERATURE_DIMENSIONS:
        raise ValueError('Invalid input units "{0}" for temperature.'.format(quantity.units))
    if outputDimensionality.simplified not in TEMPERATURE_DIMENSIONS:
        raise ValueError('Invalid output units "{0}" for temperature.'.format(units))
    
    if inputDimensionality == pq.degF.dimensionality:
        quantity = quantity + 459.67 * pq.degF
    elif inputDimensionality == pq.degC.dimensionality:
        quantity = quantity + 273.15 * pq.degC
    
    quantity = quantity.rescale(units)
    
    if outputDimensionality == pq.degF.dimensionality:
        quantity = quantity - 459.67 * pq.degF
    elif outputDimensionality == pq.degC.dimensionality:
        quantity = quantity - 273.15 * pq.degC
        
    return quantity

################################################################################

TIME_DIMENSIONS = [pq.s.simplified.dimensionality]

def convertTime(quantity, units):
    """
    Convert a given `quantity` -- a :class:`Quantity` object with units of 
    time -- to the given `units` of time. A :class:`ValueError` is raised if
    this conversion is not successful.
    """
    if isinstance(quantity, tuple):
        quantity = pq.Quantity(quantity[0], quantity[1])
    elif not isinstance(quantity, pq.Quantity):
        raise ValueError('Invalid value "{0}" for quantity; must be a tuple or Quantity object with units of time.'.format(quantity))
    
    if isinstance(units, str):
        units = pq.Quantity(1.0, units)
    
    inputDimensionality = quantity.units.dimensionality
    outputDimensionality = units.dimensionality

    if inputDimensionality.simplified not in TIME_DIMENSIONS:
        raise ValueError('Invalid input units "{0}" for time.'.format(quantity.units))
    if outputDimensionality.simplified not in TIME_DIMENSIONS:
        raise ValueError('Invalid output units "{0}" for time.'.format(units))
        
    return quantity.rescale(units)

################################################################################

LENGTH_DIMENSIONS = [pq.m.simplified.dimensionality]

def convertLength(quantity, units):
    """
    Convert a given `quantity` -- a :class:`Quantity` object with units of 
    length -- to the given `units` of length. A :class:`ValueError` is raised
    if this conversion is not successful.
    """
    if isinstance(quantity, tuple):
        quantity = pq.Quantity(quantity[0], quantity[1])
    elif not isinstance(quantity, pq.Quantity):
        raise ValueError('Invalid value "{0}" for quantity; must be a tuple or Quantity object with units of length.'.format(quantity))
    
    if isinstance(units, str):
        units = pq.Quantity(1.0, units)
    
    inputDimensionality = quantity.units.dimensionality
    outputDimensionality = units.dimensionality

    if inputDimensionality.simplified not in LENGTH_DIMENSIONS:
        raise ValueError('Invalid input units "{0}" for length.'.format(quantity.units))
    if outputDimensionality.simplified not in LENGTH_DIMENSIONS:
        raise ValueError('Invalid output units "{0}" for length.'.format(units))
        
    return quantity.rescale(units)

################################################################################

MASS_DIMENSIONS = [pq.kg.simplified.dimensionality]
MOLAR_MASS_DIMENSIONS = [(pq.kg / pq.mol).simplified.dimensionality]

def convertMass(quantity, units):
    """
    Convert a given `quantity` -- a :class:`Quantity` object with units of 
    mass -- to the given `units` of mass. A :class:`ValueError` is raised
    if this conversion is not successful. This function can handle conversion
    between intensive (molar) and extensive masses by multiplying or dividing
    by an Avogadro number as appropriate.
    """
    if isinstance(quantity, tuple):
        quantity = pq.Quantity(quantity[0], quantity[1])
    elif not isinstance(quantity, pq.Quantity):
        raise ValueError('Invalid value "{0}" for quantity; must be a tuple or Quantity object with units of mass.'.format(quantity))
    
    if isinstance(units, str):
        units = pq.Quantity(1.0, units)
        
    inputDimensionality = quantity.units.dimensionality.simplified
    outputDimensionality = units.dimensionality.simplified
    
    if inputDimensionality not in MASS_DIMENSIONS and inputDimensionality not in MOLAR_MASS_DIMENSIONS:
        raise ValueError('Invalid input units "{0}" for {1}.'.format(quantity.units, unittype))
    if outputDimensionality not in MASS_DIMENSIONS and outputDimensionality not in MOLAR_MASS_DIMENSIONS:
        raise ValueError('Invalid output units "{0}" for {1}.'.format(units, unittype))
    
    if inputDimensionality in MASS_DIMENSIONS and outputDimensionality in MOLAR_MASS_DIMENSIONS:
        quantity = quantity * constants.Na / pq.mol
    elif inputDimensionality in MOLAR_MASS_DIMENSIONS and outputDimensionality in MASS_DIMENSIONS:
        quantity = quantity / (constants.Na / pq.mol)
    
    return quantity.rescale(units)
