#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RPMDrate - Ring polymer molecular dynamics simulations
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
This module contains several simple models for interpolation of a set of
one-dimensional data, including linear, semi-logarithmic, and logarithmic
methods.
"""

import math
import numpy

################################################################################

class LinearInterpolator(object):

    def __init__(self, xdata, ydata):
        self.xdata = xdata
        self.ydata = ydata
        self.N = xdata.shape[0]
        self.xmin = self.xdata[0]
        self.xmax = self.xdata[self.N-1]
        
    def __call__(self, x):
        return self.value(x)
    
    def value(self, x):
        
        if x < self.xmin:
            raise ValueError('The requested value x = {0:g} is below the interpolation range.'.format(x))
        elif x > self.xmax:
            raise ValueError('The requested value x = {0:g} is above the interpolation range.'.format(x))
            
        xdata = self.xdata
        ydata = self.ydata
        
        # Binary search for insert position
        imin = 0; imax = self.N
        while imax >= imin:
            imid = (imin + imax) / 2
            if xdata[imid] < x:
                imin = imid + 1
            elif xdata[imid] > x:
                imax = imid - 1
            else:
                break
        
        if imin > imax:
            y2 = ydata[imin]
            y1 = ydata[imax]
            x2 = xdata[imin]
            x1 = xdata[imax]
        else:
            y2 = ydata[imax]
            y1 = ydata[imin]
            x2 = xdata[imax]
            x1 = xdata[imin]
        return (y2 - y1) / (x2 - x1) * (x - x1) + y1

################################################################################

class SemiLogXInterpolator(object):

    def __init__(self, xdata, ydata):
        self.xdata = numpy.log(xdata)
        self.ydata = ydata
        self.N = xdata.shape[0]
        self.xmin = self.xdata[0]
        self.xmax = self.xdata[self.N-1]
        
    def __call__(self, x):
        return self.value(x)
    
    def value(self, x):
        
        x = math.log(x)
        
        if x < self.xmin:
            raise ValueError('The requested value x = {0:g} is below the interpolation range.'.format(exp(x)))
        elif x > self.xmax:
            raise ValueError('The requested value x = {0:g} is above the interpolation range.'.format(exp(x)))
            
        xdata = self.xdata
        ydata = self.ydata
        
        # Binary search for insert position
        imin = 0; imax = self.N
        while imax >= imin:
            imid = (imin + imax) / 2
            if xdata[imid] < x:
                imin = imid + 1
            elif xdata[imid] > x:
                imax = imid - 1
            else:
                break
        
        if imin > imax:
            y2 = ydata[imin]
            y1 = ydata[imax]
            x2 = xdata[imin]
            x1 = xdata[imax]
        else:
            y2 = ydata[imax]
            y1 = ydata[imin]
            x2 = xdata[imax]
            x1 = xdata[imin]
        return (y2 - y1) / (x2 - x1) * (x - x1) + y1

#################################################################################

class SemiLogYInterpolator(object):

    def __init__(self, xdata, ydata):
        self.xdata = xdata
        self.ydata = numpy.log(ydata)
        self.N = xdata.shape[0]
        self.xmin = self.xdata[0]
        self.xmax = self.xdata[self.N-1]
        
    def __call__(self, x):
        return self.value(x)
    
    def value(self, x):
        
        if x < self.xmin:
            raise ValueError('The requested value x = {0:g} is below the interpolation range.'.format(x))
        elif x > self.xmax:
            raise ValueError('The requested value x = {0:g} is above the interpolation range.'.format(x))
            
        xdata = self.xdata
        ydata = self.ydata
        
        # Binary search for insert position
        imin = 0; imax = self.N
        while imax >= imin:
            imid = (imin + imax) / 2
            if xdata[imid] < x:
                imin = imid + 1
            elif xdata[imid] > x:
                imax = imid - 1
            else:
                break
        
        if imin > imax:
            y2 = ydata[imin]
            y1 = ydata[imax]
            x2 = xdata[imin]
            x1 = xdata[imax]
        else:
            y2 = ydata[imax]
            y1 = ydata[imin]
            x2 = xdata[imax]
            x1 = xdata[imin]
        return math.exp((y2 - y1) / (x2 - x1) * (x - x1) + y1)

#################################################################################

class LogLogInterpolator(object):

    def __init__(self, xdata, ydata):
        self.xdata = numpy.log(xdata)
        self.ydata = numpy.log(ydata)
        self.N = xdata.shape[0]
        self.xmin = self.xdata[0]
        self.xmax = self.xdata[self.N-1]
        
    def __call__(self, x):
        return self.value(x)
    
    def value(self, x):
        
        x = math.log(x)
        
        if x < self.xmin:
            raise ValueError('The requested value x = {0:g} is below the interpolation range.'.format(x))
        elif x > self.xmax:
            raise ValueError('The requested value x = {0:g} is above the interpolation range.'.format(x))
            
        xdata = self.xdata
        ydata = self.ydata
        
        # Binary search for insert position
        imin = 0; imax = self.N
        while imax >= imin:
            imid = (imin + imax) / 2
            if xdata[imid] < x:
                imin = imid + 1
            elif xdata[imid] > x:
                imax = imid - 1
            else:
                break
        
        if imin > imax:
            y2 = ydata[imin]
            y1 = ydata[imax]
            x2 = xdata[imin]
            x1 = xdata[imax]
        else:
            y2 = ydata[imax]
            y1 = ydata[imin]
            x2 = xdata[imax]
            x1 = xdata[imin]
        return math.exp((y2 - y1) / (x2 - x1) * (x - x1) + y1)
