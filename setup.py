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

try:
    import numpy
    from numpy.distutils.core import setup, Extension, Command
except ImportError:
    print('The numpy package is required to install RPMDrate.')
    quit()

################################################################################

class test(Command):
    """
    Run the unit test suite. All files in the `tests` directory and its
    subdirectories that match the pattern ``*test.py`` will have their unit
    test suites executed.
    """

    description = "Run the unit test suite"
    user_options = []

    def initialize_options(self):
        self.verbosity = 2

    def finalize_options(self):
        try:
            self.verbosity = int(self.verbosity)
        except ValueError:
            raise ValueError('Verbosity must be an integer.')

    def run(self):
        import sys
        if sys.version.startswith('2.6') or sys.version.startswith('3.1'):
            import unittest2 as unittest
        else:
            import unittest
        suite = unittest.TestLoader().discover('tests', pattern='*test.py', top_level_dir='.')
        unittest.TextTestRunner(verbosity=self.verbosity).run(suite)

################################################################################

# The Fortran extension modules to build using f2py
ext_modules = [
    Extension('rpmdrate._surface', ['rpmdrate/_surface.f90']),
    Extension('rpmdrate._main', ['rpmdrate/_main.pyf', 'rpmdrate/_math.f90', 'rpmdrate/_surface.f90', 'rpmdrate/_main.f90', 'rpmdrate/blas_lapack.f90'], libraries=['fftw3']),
]

setup(
    name = 'RPMDrate',
    version = '0.1.0',
    description = 'Ring polymer molecular dynamics simulations',
    author = 'Joshua W. Allen, Yury V. Suleimanov, William H. Green',
    author_email = 'rpmd_dev@mit.edu',
    url = 'http://github.com/GreenGroup/RPMDrate',
    packages = ['rpmdrate'],
    cmdclass = {'test': test},
    ext_modules = ext_modules,
    requires = ['numpy (>=1.5.0)'],
    provides = ['rpmdrate'],
)
