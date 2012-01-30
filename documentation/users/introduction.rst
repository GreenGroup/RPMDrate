********************
Introduction to RPMD
********************

What is RPMD?
=============

RPMD is short for *ring polymer molecular dynamics*: a sampling method used to
compute the approximate values of chemical dynamics properties, such as
chemical reaction rate coefficients and diffusion coefficients, that proceed
over a potential barrier. It is an extension of path integral molecular 
dynamics, which is used to compute the exact values of static chemical 
properties. RPMD exploits the isomorphism between the quantum statistical 
mechanics of a single particle and the classical statistical mechanics of a 
fictitious ring polymer composed of many copies of the original system 
connected by harmonic springs. In effect, RPMD is classical molecular dynamics
in an extended phase space, which gives the method several desirable features:

* The method becomes exact in the high temperature (classical) limit, where the
  ring polymer collapses to a single bead.

* The method gives the exact quantum mechanical result for a parabolic barrier
  bilinearly coupled to a bath of harmonic oscillators, including for 
  temperatures in the deep quantum tunneling regime.

* The method has a well-defined short-time limit that provides an upper bound
  on the classical result.

* The computed properties are rigorously independent of the choice of the
  transition state.

In short, the method captures both the quantum mechanical effects of
tunneling through a potential barrier and zero-point energy.
  
RPMD is also the name of this software package, which provides functionality
for conducting ring polymer molecular dynamics simulations in order to compute
the chemical dynamics properties of a wide variety of systems. The RPMD 
codebase is split into a Fortran 95 core used to efficiently conduct the RPMD
simulations, and a Python layer that provides a more user-friendly, scriptable
interface. This means that we don't have to choose between the flexibility of
Python and the efficiency of Fortran -- we get the best of both worlds!



Features
========



License
=======

RPMD is distributed under the terms of the 
`MIT license <http://www.opensource.org/licenses/mit-license>`_. This is a
permissive license, meaning that you are generally free to use, redistribute,
and modify RPMD at your discretion. If you use all or part of RPMD in your
work, we only ask that proper attribution be given (in addition to following
the terms of the license as stated below).

The MIT license is reproduced in its entirety below::

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

