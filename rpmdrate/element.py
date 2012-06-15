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
This module contains reference information about the chemical elements. In
particular, the atomic masses of common isotopes of many elements are
available in the ``atomicMass`` dict.
"""

################################################################################

# The atomic masses in amu of the common isotopes of many elements
atomicMass = {
    "Mu":    0.113,
    "H":     1.00782503207,
    "D":     2.0141017778,
    "T":     3.0160492777,
    "He3":   3.0160293191,
    "He":    4.00260325415,
    "Li6":   6.015122795,
    "Li":    7.01600455,
    "Be":    9.0121822,
    "B5":    10.0129370,
    "B":     11.0093054,
    "C":     12.0,
    "C13":   13.0033548378,
    "N":     14.0030740048,
    "N15":   15.0001088982,
    "O":     15.99491461956,
    "O18":   17.9991610,
    "F":     18.99840322,
    "Ne":    19.9924401754,
    "Ne22":  21.991385114,
    "Na":    22.9897692809,
    "Mg":    23.985041700,
    "Mg25":  24.98583692,
    "Mg26":  25.982592929,
    "Al":    26.98153863,
    "Si":    27.9769265325,
    "Si29":  28.976494700,
    "Si30":  29.97377017,
    "P":     30.97376163,
    "S":     31.97207100,
    "S34":   33.96786690,
    "Cl":    34.96885268,
    "Cl37":  36.96590259,
    "Ar":    39.9623831225,
    "K":     38.96370668,
    "K41":   40.96182576,
    "Ca":    39.96259098,
    "Sc":    44.9559119,
    "Ti":    47.9479463,
    "V":     50.9439595,
    "Cr50":  49.9460442,
    "Cr":    51.9405075,
    "Cr53":  52.9406494,
    "Cr54":  53.9388804,
    "Mn":    54.9380451,
    "Fe54":  53.9396105,
    "Fe":    55.9349375,
    "Fe57":  56.9353940,
    "Co":    58.9331950,
    "Ni":    57.9353429,
    "Cu":    62.9295975,
    "Cu65":  64.9277895,
    "Zn64":  63.9291422,
    "Zn66":  65.9260334,
    "Zn67":  66.9271273,
    "Zn68":  67.9248442,
    "Ga":    68.9255736,
    "Ga71":  70.9247013,
    "Ge70":  69.9242474,       
    "Ge72":  71.9220758,
    "Ge73":  72.9234589,
    "Ge74":  73.9211778,
    "Ge76":  75.9214026,
    "As":    74.9215965,
    "Se74":  73.9224764,
    "Se76":  75.9192136,
    "Se77":  76.9199140,
    "Se78":  77.9173091,
    "Se80":  79.9165213,
    "Se82":  81.9166994,
    "Br":    78.9183371,
    "Br81":  80.9162906,
    "Kr":    83.911507,
    "Rb":    84.911789738,
    "Rb87":  86.909180527,
    "Sr":    87.9056121,
    "Y":     88.9058483,
    "Zr":    89.9047044,
    "Zr91":  90.9056458,
    "Zr92":  91.9050408,
    "Zr94":  93.9063152,
    "Zr96":  95.9082734,
}
