.. _external_potenial:

****************************************
Exposing an external potential to Python
****************************************

One of the larger computational expenses in an RPMD calculation is the
evaluation of the external potential for a system. Most of the time, this
potential takes the form of a function in a compiled language, e.g. Fortran.
In this section, we demonstrate how to expose a global potential function
defined as a Fortran code to Python, and, by extension, RPMDrate.

Define the Fortran interface
============================

.. highlight:: fortran

The first step is to define the functions that Python will invoke. The first is
the function ``get_potential()`` that computes the potential and forces for
a given configuration. The parameters for ``get_potential()`` are as follows:

* ``q`` - The :math:`3 \times N_\mathrm{atoms} \times N_\mathrm{beads}` array
  of coordinates for each atom in each bead in atomic units.

* ``Natoms`` - The number of atoms in the system.

* ``Nbeads`` - The number of beads in the system.

* ``V`` - The computed potential for each bead in atomic units.

* ``dVdq`` - The :math:`3 \times N_\mathrm{atoms} \times N_\mathrm{beads}` 
  array of computed forces for each atom in each bead in atomic units.

* ``info`` - An output flag; use a nonzero value to indicate that the
  coordinates are outside of the valid range of the potential function. This
  will generally cause the RPMD trajectory to be restarted using the last
  previous known good coordinates.

A template for ``get_potential()`` is given below::

    subroutine get_potential(q, Natoms, Nbeads, V, dVdq, info)
    
        implicit none
        integer, intent(in) :: Natoms
        integer, intent(in) :: Nbeads
        double precision, intent(in) :: q(3,Natoms,Nbeads)
        double precision, intent(out) :: V(Nbeads)
        double precision, intent(out) :: dVdq(Nbeads)
        integer, intent(out) :: info
        
        integer :: k
        
        ! Iterate over the beads, computing the potential and forces of each separately
        do k = 1, Nbeads
        
            ! Fill in the code for calling the Fortran potential function here.
            ! Use q(:,:,k) for the position
            ! Place the computed potential in V(k) and forces in dVdq(:,:,k) 
        
        end do
        
    end subroutine

Many potential energy functions require some initialization. If yours is one
of these, then you must write a second function ``initialize_potential()`` that
RPMDrate can call to perform the initialization. A template for this function 
is given below::

    subroutine initialize_potential()
    
        ! Fill in the code for initializing the Fortran potential function here.
        
    end subroutine

Define the f2py wrapper module
==============================

The Fortran functions defined above are exposed to Python using the 
`f2py <http://www.scipy.org/F2py>`_ utility. This utility has been shipped with
NumPy since 2007, and should therefore not require any additional setup if you
are using a recent version of NumPy. Note that you will need a reasonably
recent Fortran compiler, such as `gfortran <http://gcc.gnu.org/fortran/>`_, 
`g95 <http://www.g95.org/>`_, or 
`ifort <http://software.intel.com/en-us/articles/intel-compilers/>`_.

If you have created the ``get_potential()`` and ``initialize_potential()``
Fortran subroutines as shown in the previous section, then the following f2py
module definition file should suffice::

    python module PES
        interface
            subroutine initialize_potential()
            end subroutine initialize_potential
            subroutine get_potential(q,Natoms,Nbeads,V,dVdq,info)
                integer, intent(in) :: Natoms
                integer, intent(in) :: Nbeads
                double precision dimension(1:3,1:Natoms,1:Nbeads), intent(in) :: q
                double precision dimension(1:Nbeads), intent(out) :: V
                double precision dimension(1:3,1:Natoms,1:Nbeads), intent(out) :: dVdq
                integer, intent(out) :: info
            end subroutine get_potential
        end interface 
    end python module PES

Place the above in a file called ``PES.pyf`` in the same directory as the
RPMDrate input file and the other Fortran source files for your potential.

Compile the Python wrapper
==========================

.. highlight:: bash

All that remains is to perform the actual compilation. For the purposes of this
example, we assume that the entirely of the external potential is divided into
two files:

* ``PES.f`` - The Fortran source file containing the potential, without any
  modifications for use by RPMDrate.

* ``PES_wrapper.f90`` - The Fortran source file containing the 
  ``get_potential()`` and ``initialize_potential()`` Fortran subroutines that
  interface between RPMDrate and the potential function.

With these two files and the f2py module definition file ``PES.pyf``, we can
invoke the compilation using the following command::

    $ f2py -c PES.f PES_wrapper.f90 PES.pyf

.. highlight:: python

You will see a large amount of compiler messages scroll by, including possibly
some warnings depending on how compliant your Fortran potential function is.
If the compilation is successful, you should see a new file -- ``PES.so`` on
Mac/Linux/etc., ``PES.pyd`` on Windows -- that represents the wrapped Fortran
code. This file is ready for use by RPMDrate; this line from the input file
tells RPMDrate about the ``get_potential()`` and ``initialize_potential()``
subroutines::  

    from PES import get_potential, initialize_potential
