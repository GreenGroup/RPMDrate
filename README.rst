**************************************************
RPMD - Ring polymer molecular dynamics simulations
**************************************************

About RPMD
==========

RPMD is a free, open-source software package for using ring polymer molecular
dynamics simulations to compute properties of chemical systems, including, but
not limited to, chemical reaction rates.

License
=======

RPMD is distributed under the terms of the 
`MIT license <http://www.opensource.org/licenses/mit-license>`_::

    Copyright (c) 2012 by Joshua W. Allen (jwallen@mit.edu)
                          Yury V. Suleimanov (ysuleyma@mit.edu)
                          William H. Green (whgreen@mit.edu)
    
    Permission is hereby granted, free of charge, to any person obtaining a 
    copy of this software and associated documentation files (the "Software"), 
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense, 
    and/or sell copies of the Software, and to permit persons to whom the 
    Software is furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
    DEALINGS IN THE SOFTWARE. 

Installation
============

Dependencies
------------

RPMD depends on several other packages in order to provide its full
functional capabilities. The following dependencies are required to compile
and run RPMD:

* `Python <http://www.python.org/>`_ (version 2.5.x or later, including any version of Python 3, is recommended)

* `NumPy <http://numpy.scipy.org/>`_ (version 1.5.0 or later is recommended)

* `FFTW <http://www.fftw.org/>`_ (version 3.3 or later is recommended)

* A standard Fortran 90/95 compiler (e.g. ``gfortran``, ``g95``, ``ifort``, etc.)

RPMD uses the ``f2py`` tool (bundled with NumPy) to expose its Fortran code to
Python.

Getting RPMD
------------

The best way to obtain a copy of the repository is to clone it using `git
<http://git-scm.com/>`_::

    $ git clone git://github.com/GreenGroup/RPMD.git

This enables you to easy update your local clone with the latest changes. If
you intend to contribute, you should fork the project on GitHub first, and
clone from your fork.

Compiling from Source
---------------------

RPMD is installed by invoking the ``setup.py`` script::

    $ python setup.py install

This will compile the Fortran source and Python wrapper code, which may take
some time. Note that you may need administrator privileges to install RPMD.

If you wish to use RPMD without installing, simply add the folder containing
this file to your ``PYTHONPATH`` environment variable and compile the source
code in-place using the command ::

    $ python setup.py build_ext --inplace

A Makefile that wraps these commands has been provided. The Makefile also
provides a clean target for deleting all files created during compilation.

Running the Unit Tests
----------------------

RPMD comes with a large suite of unit tests that ensure functionality is
working as intended. To run these tests, first install RPMD or compile it from
source in-place using the directions above. Then, simply invoke the entire
suite of unit tests using ::

    $ python setup.py test

This will run all of the unit tests in sequence, which may take some time. If
one or more unit tests fail, please report them on the GitHub issue tracker
(if they aren't already reported).
