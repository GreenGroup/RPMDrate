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

r"""
This module contains classes and methods for working with physical constants
relevant in chemistry applications. The available physical constants are listed
in the table below:

.. table:: Physical constants defined in the :mod:`rpmd.constants` module
    
    =========== ======================= =========================================================== ================================================================
    Constant    Symbol                  Value                                                       Description
    =========== ======================= =========================================================== ================================================================
    ``Na``      :math:`N_\mathrm{A}`    :math:`6.02214179 \times 10^{23} \ \mathrm{mol^{-1}}`       Avogadro constant
    ``R``       :math:`R`               :math:`8.314472 \ \mathrm{J/mol \cdot K}`                   gas law constant
    ``c``       :math:`c`               :math:`299792458 \ \mathrm{m/s}`                            speed of light in a vacuum
    ``e``       :math:`e`               :math:`1.602176565 \times 10^{-19} \ \mathrm{C}`            elementary charge
    ``h``       :math:`h`               :math:`6.62606896 \times 10^{-34} \ \mathrm{J \cdot s}`     Planck constant
    ``hbar``    :math:`\hbar`           :math:`1.054571726 \times 10^{-34} \ \mathrm{J \cdot s}`    reduced Planck constant
    ``kB``      :math:`k_\mathrm{B}`    :math:`1.3806504 \times 10^{-23} \ \mathrm{J/K}`            Boltzmann constant
    ``pi``      :math:`\pi`             :math:`3.14159 \ldots`      
    =========== ======================= =========================================================== ================================================================

The physical constants are defined as module-level variables. Accordingly, the
recommended method of importing this module is ::

    import rpmd.constants as constants

so as to not place the constants in the importing module's global namespace.
"""

import math

################################################################################

#: :math:`\pi = 3.14159 \ldots`      
pi = float(math.pi)

#: The Avogadro constant :math:`N_\mathrm{A}` in :math:`\mathrm{mol^{-1}}`       
Na = 6.02214179e23

#: The Boltzmann constant :math:`k_\mathrm{B}` in :math:`\mathrm{J/K}`
kB = 1.3806504e-23

#: The gas law constant :math:`R` in :math:`\mathrm{J/mol \cdot K}`
R = 8.314472

#: The Planck constant :math:`h` in :math:`\mathrm{J \cdot s}`
h = 6.62606896e-34

#: The reduced Planck constant :math:`\hbar` in :math:`\mathrm{J \cdot s}`
hbar = 1.054571726e-34

#: The speed of light in a vacuum :math:`c` in :math:`\mathrm{m/s}`
c = 299792458

#: The elementary charge :math:`e = 1.602176565 \times 10^{-19}` in :math:`\mathrm{C}`
e = 1.602176565e-19
