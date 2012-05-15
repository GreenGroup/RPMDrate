******************************
Browsing RPMDrate output files
******************************

RPMDrate generates a large number of output files for each calculation. These
are sorted into subfolders by temperature and number of beads, so that the
same input file can be used by multiple RPMD calculations. RPMDrate will also
use these output files for restarting incomplete jobs, to avoid repeating
already-finished steps in the calculation.

Umbrella configurations
=======================

The computed umbrella configurations are saved to
``umbrella_configurations.dat`` in the top-level directory, as the 
configurations are independent of the temperature and number of beads.
This file contains the Cartesian coordinates in atomic units for each
configuration.

Umbrella sampling
=================

The computed values for the mean :math:`\left< \xi \right>` and variance
:math:`\sigma_\xi^2` of the reaction coordinate are saved to a set of files
named ``umbrella_sampling_*.dat``, where ``*`` represents the value of the
reaction coordinate :math:`\xi`.

Potential of mean force
=======================

The potential of mean force as a function of reaction coordinate, as computed
by umbrella integration, is saved to the file ``potential_of_mean_force.dat``.

Recrossing factor
=================

The computed value of the recrossing factor is saved to a file named
``recrossing_factor_*.dat``, where ``*`` represents the value of the reaction
coordinate :math:`\xi` at which the recrossing factor was computed.

Rate coefficient
================

The final value of the rate coefficient, along with a summary of the various
intermediate values used to compute the rate coefficient, is saved to a file
named ``rate_coefficient_*.dat``, where ``*`` represents the value of the
reaction coordinate :math:`\xi` at which the recrossing factor was computed.
Note that, for a converged calculation, the same value of the rate coefficient
should be obtained using any value of :math:`\xi` for computing the potential
of mean force and recrossing factor.
 