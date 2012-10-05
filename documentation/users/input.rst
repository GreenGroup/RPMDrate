*******************************
Creating an RPMDrate input file
*******************************

An RPMDrate input file is simply a Python script.

Define the potential function
=============================

The RPMDrate code must be provided with an external function that computes the
potential for a given geometry. In principle, this can be any Python callable
object of the form ::

    V, dVdq = get_potential(q)

where the lone parameter ``q`` is the :math:`3 \times N_\mathrm{atoms} \times N_\mathrm{beads}`
array of positions in atomic units, and the outputs are the corresponding
potentials and forces for each bead. The input file requires the name to be
``get_potential``.

In practice, evaluating the potential and forces is often computationally
expensive, and so the potential function is likely written in another language
(e.g. Fortran) and must therefore be exposed to Python. The process of doing
this is explained in :ref:`external_potenial`.

If your implementation of the potential function requires any initialization,
you can also provide a Python function ``initialize_potential`` for this
purpose. This function is assumed to not accept any parameters, and any return
value is ignored. If the ``initialize_potential`` identifier is not found,
the code will assume that no initialization of the potential is required.

A recommended practice is to define the potential function (and its
initialization subroutine, if applicable) in a separate module ``PES`` placed
in the same directory as the input file. If this is done correctly, then the
input file need only import them::

    from PES import get_potential, initialize_potential

If no initialization subroutine is required, then the syntax is :: 

    from PES import get_potential

Define the reactants
====================

The next step is to provide information about the bimolecular reactants. This
is done by calling the ``reactants()`` function with the following arguments,
all of which are required:

* ``atoms`` - A list of the atom symbols for all atoms in the reactant pair.
  RPMDrate recognizes a variety of atom symbols; the full list is available in
  the :mod:`rpmd.element` module.

* ``reactant1Atoms`` - A list of the indices of each atom in the first reactant
  molecule. The indices correspond to those of the ``atoms`` list, except that
  the first index is one, not zero.

* ``reactant2Atoms`` - A list of the indices of each atom in the second reactant
  molecule. The indices correspond to those of the ``atoms`` list, except that
  the first index is one, not zero.

* ``Rinf`` - The center-of-mass distance at which interaction between the
  reactant molecules is negligible. This is specified using a 2-tuple of the
  form ``(value,units)``, where ``value`` is a number and ``units`` is a
  string with the corresponding units (of length). RPMDrate can handle a
  variety of input units for length, including ``m``, ``angstrom``, and 
  ``bohr`` (atomic units).

The following example describes the bimolecular reaction :math:`\mathrm{CH_4 + H} \rightarrow \mathrm{CH_3 + H_2}`::
    
    reactants(
        atoms = ['C', 'H', 'H', 'H', 'H', 'H'],
        reactant1Atoms = [1,2,3,4,5],
        reactant2Atoms = [6],
        Rinf = (15,"angstrom"),
    )

Define the transition state
===========================

Next we define parameters for RPMDrate to use to determine the transition
state dividing surface. The easiest way to do this is to specify a transition
state geometry and lists of the bonds that are forming and breaking. RPMDrate
will use this information to extract the relevant transition state bond lengths
to use in the transition state dividing surface. 

* The ``geometry`` is specified as a 2-tuple of the form ``(value,units)``, 
  where ``value`` is a :math:`N_\mathrm{atoms} \times 3` matrix (list of lists)
  of coordinates and ``units`` is a string with the corresponding units (of
  length). RPMDrate can handle a variety of input units for length, including
  ``m``, ``angstrom``, and ``bohr`` (atomic units). The coordinates must be
  given in the same order as the list of atom symbols ``atoms`` was given in
  the ``reactants()`` block.

* The ``formingBonds`` and ``breakingBonds`` are lists of 2-tuples containing
  the indices of the atoms involved in each forming or breaking bond. As before,
  the indices correspond to the list of atom symbols ``atoms`` was given in
  the ``reactants()`` block, and the indices start from one instead of zero.

For our :math:`\mathrm{CH_4 + H} \rightarrow \mathrm{CH_3 + H_2}` example, we
define the transition state as follows::

    transitionState(
        geometry = (
            [[  -0.00480000,  0.01160000,  1.09079000 ],
             [   1.15021585,  0.04354511,  1.61546028 ],
             [  -0.02126938,  0.04768982, -0.17758979 ],
             [  -0.81242986, -0.81401508,  1.61657965 ],
             [  -0.42859126,  0.98191535,  1.36992096 ],
             [  -0.80323584,  1.83970265,  1.61663210 ]],
            "angstrom",
        ),
        formingBonds = [(5,6)], 
        breakingBonds = [(1,5)],
    )

Some systems will have multiple equivalent reaction pathways; in our example,
the H radical can abstract any of the H atoms from methane. These equivalent
transition states can be added with one or more ``equivalentTransitionState()``
blocks. In these blocks you need only specify the forming and breaking bonds
in an order corresponding to that of the main ``transitionState()`` block.
RPMDrate will automatically copy the relevant transition state distances from
the main transition state to use in each of the equivalent transition states.

In our example we must add three equivalent transition states to obtain the
correct reaction-path degeneracy of four::
 
    equivalentTransitionState(
        formingBonds=[(4,6)], 
        breakingBonds=[(1,4)],
    )
    equivalentTransitionState(
        formingBonds=[(3,6)], 
        breakingBonds=[(1,3)],
    )
    equivalentTransitionState(
        formingBonds=[(2,6)], 
        breakingBonds=[(1,2)],
    )

Define the thermostat
=====================

In molecular dynamics simulations, a thermostat is applied to control the
system temperature by simulating a coupling to a heat bath. RPMDrate provides
a choice of two thermostats: Andersen and GLE.

Andersen thermostat
-------------------

The Andersen thermostat is a simple thermostat that periodically resamples
the momenta from a Boltzmann distribution at the desired temperature to
simulate a collision with the heat bath. The following sets up an Andersen
thermostat with an automatically-determined sampling time::

    thermostat('Andersen')

You can also set the sampling time explicitly using the ``samplingTime``
parameter. This parameter is given as a 2-tuple of the form ``(value,units)``, 
where ``value`` is a number and ``units`` is a string with the corresponding
units (of time). RPMDrate can handle a variety of input units for time, 
including ``s``, ``ms``, ``us``, ``ns``, ``ps``, and ``fs``. An example of
specifying the sampling time for the Andersen thermostat is given below::

    thermostat('Andersen', samplingTime=(0.1,'ps'))

GLE thermostat
--------------

The colored-noise generalized Langevin equation (GLE) thermostat is a more
advanced thermostat, requiring more effort to set up but also generally giving
faster convergence to the quantum mechanical result. To use this thermostat,
you must provide at least the drift matrix :math:`A` in units of 
:math:`\mathrm{s^{-1}}` in a separate text file. For example, the following
sets up a GLE thermostat with the drift matrix read from the file ``gle_A.txt``
in the same directory as the input file::

    thermostat('GLE', A='gle_A.txt')

You can optionally also provide a diffusion matrix :math:`C` in units of
:math:`\mathrm{K}` in another input file. In this case your input would look
similar to ::

    thermostat('GLE', A='gle_A.txt', A='gle_C_1000.txt')

Note that the diffusion matrix is a function of temperature, so you will need
to provide a different one at each temperature. If you do not provide a
diffusion matrix, a default version will be used instead.

The files containing the :math:`A` and :math:`C` matrices should be structured
such that each line contains one row of the matrix, with each element in that
row separated by whitespace. The hash character ``#`` can be used to indicate
comment lines. Several of the provided examples come with drift and diffusion
matrices. 

An easy way to generate the :math:`A` and :math:`C` matrices is to use the
online form at the GLE4MD web page at 
http://gle4md.berlios.de/compose.php?page=matrix. On this page, select 
``PI+GLE`` as the GLE type, and then fill out the rest of the form with the
values for your system. Be sure to check that the output units for the two
matrices are :math:`\mathrm{s^{-1}}` and :math:`\mathrm{K}`, respectively. 

Generate the umbrella configurations
====================================

At this point we have successfully set up the system parameters for the RPMD
rate coefficient calculation. In the remainder of the input file we tell
RPMDrate the parameters to use for the various steps involved in computing the
rate coefficient for the input system.

The first step is to generate the initial configurations required for umbrella
sampling and integration. This is accomplished by a brief umbrella sampling
in each window using only one bead, using the result from the previous window
as the initial configuration for the next window. The resulting configurations
are independent of temperature and number of beads, and can therefore be
reused in subsequent calculations.

The generation of the umbrella configurations is invoked using a
``generateUmbrellaConfigurations()`` block, which accepts the following 
parameters:

* ``dt`` - The time step size, as a 2-tuple of the form ``(value,units)``, 
  where ``value`` is a number and ``units`` the corresponding units (of time).
  RPMDrate recognizes many units of time, including s, ms, us, ns, ps, and fs.

* ``evolutionTime`` - The total length of each trajectory, as a 2-tuple of the
  form ``(value,units)``, where ``value`` is a number and ``units`` the
  corresponding units (of time).

* ``xi_list`` - A list, tuple, or array containing the values of the reaction
  coordinate to generate the umbrella configurations at.

* ``kforce`` - The value of the umbrella integration force constant to use,
  in atomic units.

An example of a ``generateUmbrellaConfigurations()`` block is given below::

    generateUmbrellaConfigurations(
        dt = (0.0001,"ps"),
        evolutionTime = (5,"ps"),
        xi_list = numpy.arange(-0.05, 1.15, 0.01),
        kforce = 0.1 * T,
    )

Define the umbrella integration windows
=======================================

The next step is to define the umbrella integration windows. As the RPMDrate
input file is a Python script, this is very easily done with a small bit of
Python code. The idea is to create a list of the umbrella integration windows,
each one a ``Window()`` object with the following parameters: 

* ``xi`` - The value of the reaction coordinate for the window.

* ``kforce`` - The value of the umbrella integration force constant to use for
  the window, in atomic units.

* ``trajectories`` - The number of distinct umbrella sampling trajectories to
  run for the window. Each of these trajectores can be run in parallel.
  
* ``equilibrationTime`` - The amount of simulation time spent equilibrating
  each trajectory before sampling is initiated. Equilibration allows for
  randomization of the initial position and momentum, which gives better
  statistics in the sampling. This value is specified as a 2-tuple of the
  form ``(value,units)``, where ``value`` is a number and ``units`` the
  corresponding units (of time).

* ``evolutionTime`` - The total sampling length of each trajectory, as a 
  2-tuple of the form ``(value,units)``, where ``value`` is a number and
  ``units`` the corresponding units (of time).

An example of code that generates a set of umbrella integration windows is
given below::

    windows = []
    for xi in numpy.arange(-0.05, 1.15, 0.01):
        window = Window(xi=xi, kforce=0.1*T, trajectories=100, equilibrationTime=(20,"ps"), evolutionTime=(100,"ps"))
        windows.append(window)

Note that you do not necessarily have to space out the windows evenly, or to
use the same force constant or sampling times in each window. For example,
the following code uses a larger reaction coordinate step size at lower values
of ``xi``::

    windows = []
    for xi in numpy.arange(-0.05, 0.25, 0.05):
        window = Window(xi=xi, kforce=0.1*T, trajectories=100, equilibrationTime=(20,"ps"), evolutionTime=(100,"ps"))
        windows.append(window)
    for xi in numpy.arange(-0.25, 0.55, 0.03):
        window = Window(xi=xi, kforce=0.1*T, trajectories=100, equilibrationTime=(20,"ps"), evolutionTime=(100,"ps"))
        windows.append(window)
    for xi in numpy.arange( 0.55, 1.15, 0.01):
        window = Window(xi=xi, kforce=0.1*T, trajectories=100, equilibrationTime=(20,"ps"), evolutionTime=(100,"ps"))
        windows.append(window)

Conduct the umbrella sampling
=============================

Once we have defined the umbrella sampling windows, we can conduct the sampling
itself via a ``conductUmbrellaSampling()`` block, which accepts the following 
parameters:

* ``dt`` - The time step size, as a 2-tuple of the form ``(value,units)``, 
  where ``value`` is a number and ``units`` the corresponding units (of time).
  RPMDrate recognizes many units of time, including s, ms, us, ns, ps, and fs.

* ``windows`` - The list of umbrella sampling windows, as generated in the
  previous section.

An example of a ``conductUmbrellaSampling()`` block is given below::

    conductUmbrellaSampling(
        dt = (0.0001,"ps"),
        windows = windows,
    )
    
Compute the potential of mean force
===================================

When umbrella sampling is complete, we are ready to compute the potential of
mean force using the technique of umbrella integration. This is done using a
``computePotentialOfMeanForce()`` block, which accepts the following 
parameters:

* ``windows`` - The list of umbrella sampling windows used to conduct the
  umbrella sampling.

* ``xi_min`` - The value of the reaction coordinate to use as a lower bound
  for the umbrella integration range.

* ``xi_max`` - The value of the reaction coordinate to use as an upper bound
  for the umbrella integration range.

* ``bins`` - The number of bins for the umbrella integration algorithm to use.

An example of a ``computePotentialOfMeanForce()`` block is given below::

    computePotentialOfMeanForce(windows=windows, xi_min=0.0, xi_max=1.1, bins=5000)

Compute the recrossing factor
=============================

The final step before calculating the rate coefficient is to compute the
recrossing factor. It is the combination of this factor and the potential of
mean force that makes the RPMD technique for computing the rate coefficient
independent of the choice of dividing surface. The approach is to start a large
number of trajectories from the transition state geometry, and determine which
fraction of these trajectories result in products. A parent trajectory is
evolved in time while constrained to the transition state; occasionally a
number of short, unconstrained child trajectories are spawned to determine the
fraction that result in products.

The computation of the recrossing factor is invoked using a
``computeRecrossingFactor()`` block, which accepts the following 
parameters:

* ``dt`` - The time step size, as a 2-tuple of the form ``(value,units)``, 
  where ``value`` is a number and ``units`` the corresponding units (of time).

* ``equilibrationTime`` - The amount of simulation time spent equilibrating
  each trajectory before sampling is initiated. Equilibration allows for
  randomization of the initial position and momentum, which gives better
  statistics in the sampling. This value is specified as a 2-tuple of the
  form ``(value,units)``, where ``value`` is a number and ``units`` the
  corresponding units (of time).

* ``childTrajectories`` - The total number of child trajectories to run,
  including many samplings from the parent trajectory.

* ``childSamplingTime`` - The amount of time to evolve the parent trajectory
  between each set of child trajectories, as a 2-tuple of the form 
  ``(value,units)``, where ``value`` is a number and ``units`` the
  corresponding units (of time).

* ``childrenPerSampling`` - The number of child trajectores to sample after
  each evolution of the parent trajectory.

* ``childEvolutionTime`` - The total length of each child trajectory, as a
  2-tuple of the form ``(value,units)``, where ``value`` is a number and
  ``units`` the corresponding units (of time).

An example of a ``computeRecrossingFactor()`` block is given below::

    computeRecrossingFactor(
        dt = (0.0001,"ps"),
        equilibrationTime = (20,"ps"),
        childTrajectories = 100000,
        childSamplingTime = (2,"ps"),
        childrenPerSampling = 100,
        childEvolutionTime = (0.05,"ps"),
    )

The ``computeRecrossingFactor()`` block will automatically use the maximum of
the potential of mean force for its calculation, as this will give the largest
value of the recrossing factor. You can override this by specifying an
additional parameter ``xi_current`` with the value of the reaction coordinate
you wish to use for the recrossing factor calculation. An example of this is
given below::

    computeRecrossingFactor(
        dt = (0.0001,"ps"),
        equilibrationTime = (20,"ps"),
        childTrajectories = 100000,
        childSamplingTime = (2,"ps"),
        childrenPerSampling = 100,
        childEvolutionTime = (0.05,"ps"),
        xi_current = 1.016,
    )

Compute the rate coefficient
============================

We are now ready to compute the bimolecular rate coefficient for this system.
This is done using a ``computeRateCoefficient()`` block, as shown below::

    computeRateCoefficient()

The ``computeRateCoefficient()`` block does not require any parameters. It
uses the most recently computed potential of mean force and recrossing factor
to compute the value of the bimolecular rate coefficient. If you have only
done the potential of mean force calculation, you can still use this to
generate an estimate of the rate coefficient; in this case, a value of unity
will be assumed for the recrossing factor.

Additional input file options
=============================

This section discusses additional items you can optionally add to the input
file to modify how the RPMD calculation is performed. In all cases these items
are entirely optional.

Define the random number seed
-----------------------------

If you wish, you can specify a particular integer seed value to use to 
initialize the random number generator. This is entirely optional, as RPMDrate
will automatically seed the random number generator (based on the system clock)
unless this value is provided. To use a particular value, add a line to your
input file of the form ::

    random_seed = 1

where you replace the value ``1`` with your desired integer seed value.
