****************
Running RPMDrate
****************

.. highlight:: bash

Once you have set up the input file and exposed the Fortran potential function
to Python, you are ready to run the RPMD calculation. To do this, you invoke
the main ``rpmdrate.py`` script, passing the input file, temperature, and
number of beads as command-line arguments::

    $ python rpmdrate.py examples/H+CH4/input.py 1000 1

Specifying the number of processors
-----------------------------------

RPMDrate is also designed to exploit multiprocessor systems. To specify the
number of processors, use the ``-p`` flag. For example, the following indicates
that eight processors are available for RPMDrate to use::

    $ python rpmdrate.py examples/H+CH4/input.py 1000 1 -p 8
