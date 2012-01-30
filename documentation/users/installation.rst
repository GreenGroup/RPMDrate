***************
Installing RPMD
***************

There are multiple ways to install RPMD on your computer:

* Use one of the automatic Python package installation tools, such as ``pip``
  (recommended) or the older ``easy_install``.

* Download and compile a stable release package as a zip file or tarball.

* Clone our git repository and easily follow along with the bleeding-edge 
  development as it happens.



Using a Package Installer
=========================

The easiest way to install RPMD is using an automatic package installer tool
such as ``pip`` (recommended) or ``easy_install``::

$ pip install rpmd

or ::

$ easy_install rpmd

The benefit of this approach is that all of the dependencies will be handled
automatically for you. All you need is a reasonably recent Fortran 90/95
compiler for compiling the Fortran modules.



Compiling from Source
=====================

If you prefer, you can also download a RPMD release package manually and
install that way. You might also wish to follow the bleeding-edge development
package.

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

Compiliation
------------

Once you have obtained a copy of RPMD, either from a released package or
from cloning the git repository, RPMD can be installed by invoking the 
``setup.py`` script as usual::

$ python setup.py install

This will compile the Fortran source and Python wrapper code, which may take
some time. Note that you may need administrator privileges to install RPMD.

If you wish to use RPMD without installing, simply add the folder containing
this file to your ``PYTHONPATH`` environment variable and compile the source
code in-place using the command ::

$ python setup.py build_ext --inplace

A Makefile that wraps these commands has been provided. The Makefile also
provides a clean target for deleting all files created during compilation.



Developing for RPMD
===================

If you would like to follow along with the bleeding-edge development of RPMD,
the best way to do so is to clone the repository using our version control
system, `git <http://git-scm.com/>`_::

$ git clone git://github.com/GreenGroup/RPMD.git

This enables you to easy update your local copy with the latest changes. If
you intend to contribute to RPMD - and you are welcome to do so! - you should 
fork the project on GitHub first, and clone from your fork.

Once you have cloned the repository, follow the compilation instructions
for building RPMD in-place. You may also wish to add the folder containing the
RPMD repository to your ``PYTHONPATH`` environment variable.

Running the Unit Tests
======================

RPMD comes with a large suite of unit tests that ensure functionality is
working as intended. To run these tests, first install RPMD or compile it
from source in-place using the directions in the previous section. Then, simply
invoke the entire suite of unit tests using ::

$ python setup.py test

This will run all of the unit tests in sequence, which may take some time. If
one or more unit tests fail, please report them on the GitHub issue tracker
(if they aren't already reported).
