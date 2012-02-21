#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RPMD - Ring polymer molecular dynamics simulations
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

import math
import numpy
import logging
import quantities as pq
import multiprocessing

import rpmd.constants as constants

from rpmd._main import *

################################################################################

class RPMDError(Exception):
    """
    An exception raised when an error occurs during an RPMD simulation. Pass a
    string describing the circumstances of the exceptional behavior.
    """
    pass

################################################################################

def runUmbrellaTrajectory(rpmd, xi_current, q, equilibrationSteps, evolutionSteps, saveTrajectory):
    """
    Run an individual umbrella integration trajectory, returning the sum of the
    first and second moments of the reaction coordinate at each time step.
    """
    rpmd.activate()
    p = rpmd.sampleMomentum()
    result = system.equilibrate(0, p, q, equilibrationSteps, xi_current, rpmd.potential, False, saveTrajectory)
    dav, dav2, result = system.umbrella_trajectory(0, p, q, evolutionSteps, xi_current, rpmd.potential, saveTrajectory)
    return dav, dav2

def runRecrossingTrajectory(rpmd, xi_current, p, q, evolutionSteps, saveTrajectory):
    """
    Run an individual recrossing factor child trajectory, returning the 
    contributions to the numerator and denominator of the recrossing factor
    from this trajectory.
    """
    rpmd.activate()
    kappa_num = numpy.zeros(evolutionSteps, order='F')
    kappa_denom = numpy.array(0.0, order='F')
    result = system.recrossing_trajectory(0, p, q, xi_current, rpmd.potential, saveTrajectory, kappa_num, kappa_denom)
    return kappa_num, kappa_denom

class RPMD:
    """
    A representation of a ring polymer molecular dynamics (RPMD) job for
    computing gas-phase chemical reaction rates. The attributes are:
    
    =========================== ================================================
    Attribute                   Description
    =========================== ================================================
    `mass`                      The mass of each atom in the molecular system
    `Natoms`                    The number of atoms in the molecular system
    `reactants`                 The dividing surface near the reactants, as a :class:`Reactants` object
    `transitionState`           The dividing surface near the transition state, as a :class:`TransitionState` object
    `potential`                 A function that computes the potential and forces for a given position
    --------------------------- ------------------------------------------------
    `beta`                      The reciprocal temperature of the RPMD simulation
    `dt`                        The time step to use in the RPMD simulation
    `Nbeads`                    The number of beads per atom in the RPMD simulation
    `xi_current`                The current value of the reaction coordinate
    `mode`                      A flag indicating the type of RPMD calculation currently underway (1 = umbrella, 2 = recrossing)
    =========================== ================================================
    
    """

    def __init__(self, reactants, transitionState, potential):
        """
        Initialize an RPMD object. The `mass` of each atom should be given in
        g/mol, while the `Rinf` value should be given in angstroms. (They will
        be converted to atomic units.)
        """
        self.mass = reactants.mass
        self.Natoms = len(self.mass)
        self.reactants = reactants
        self.transitionState = transitionState
        self.potential = potential
        
        self.beta = 0
        self.dt = 0
        self.Nbeads = 0
        self.xi_current = 0
        self.mode = 0
        
    def activate(self):
        """
        Set this object as the active RPMD system in the Fortran layer. Note
        that the dividing surface properties must be set in the corresponding
        modules in ``rpmd._main``, *not* in ``rpmd._surface``.
        """
        Natoms = self.mass.shape[0]
        system.dt = self.dt
        system.beta = self.beta
        system.kforce = self.kforce
        system.mass[0:Natoms] = self.mass
        system.mode = self.mode
        self.reactants.activate(module=reactants)
        self.transitionState.activate(module=transition_state)
    
    def computeStaticFactor(self, T, Nbeads, dt, 
                            xi_list,
                            equilibrationTime,
                            numberOfTrajectories,
                            evolutionTime,
                            kforce,
                            processes=1,
                            saveTrajectories=False):
        """
        Return the value of the static factor :math:`p^{(n)}(s_1, s_0)` as
        computed using umbrella integration.
        """
        
        # Set the parameters for the RPMD calculation
        self.beta = 4.35974417e-18 / (constants.kB * T)
        self.dt = dt / 2.418884326505e-5
        self.Nbeads = Nbeads
        Nxi = len(xi_list)
        equilibrationTime /= 2.418884326505e-5
        evolutionTime /= 2.418884326505e-5
        geometry = self.transitionState.geometry[:,:,0]
        self.kforce = kforce
        self.mode = 1
        
        av = numpy.zeros(Nxi)
        av2 = numpy.zeros(Nxi)
        
        # Create a pool of subprocesses to farm out the individual trajectories to
        pool = multiprocessing.Pool(processes=processes)
        results = []

        logging.info('******************')
        logging.info('RPMD static factor')
        logging.info('******************')
        logging.info('')
        
        equilibrationSteps = int(round(equilibrationTime / self.dt))
        evolutionSteps = int(round(evolutionTime / self.dt))
        
        logging.info('Parameters')
        logging.info('==========')
        logging.info('Temperature                             = {0:g} K'.format(T))
        logging.info('Number of beads                         = {0:d}'.format(Nbeads))
        logging.info('Time step                               = {0:g} ps'.format(self.dt * 2.418884326505e-5))
        logging.info('Number of umbrella integration windows  = {0:d}'.format(Nxi))
        logging.info('Initial equilibration time              = {0:g} ps ({1:d} steps)'.format(equilibrationSteps * self.dt * 2.418884326505e-5, equilibrationSteps))
        logging.info('Trajectory evolution time               = {0:g} ps ({1:d} steps)'.format(evolutionSteps * self.dt * 2.418884326505e-5, evolutionSteps))
        logging.info('Number of trajectories per window       = {0:d}'.format(numberOfTrajectories))
        logging.info('')

        self.activate()

        # Generate initial position using transition state geometry
        # (All beads start at same position)
        q = numpy.zeros((3,self.Natoms,self.Nbeads), order='F')
        for i in range(3):
            for j in range(self.Natoms):
                for k in range(self.Nbeads):
                    q[i,j,k] = geometry[i,j]

        # Find the window nearest to the dividing surface
        for start in range(Nxi):
            if xi_list[start] >= 1:
                break
           
        # Equilibrate in each window to determine the initial positions
        # First start at xi = 1 and move in the xi > 1 direction, using the
        # result of the previous xi as the initial position for the next xi
        q_initial = numpy.zeros((3,self.Natoms,self.Nbeads,Nxi), order='F')
        for l in range(start, Nxi):
            xi_current = xi_list[l]
            
            # Equilibrate in this window
            logging.info('Equilibrating trajectory at xi = {0:g} for {1:g} ps...'.format(xi_current, equilibrationSteps * self.dt * 2.418884326505e-5))
            p = self.sampleMomentum()
            result = system.equilibrate(0, p, q, equilibrationSteps, xi_current, self.potential, False, saveTrajectories)
            logging.info('Finished equilibrating trajectory at xi = {0:g}.'.format(xi_current))
            q_initial[:,:,:,l] = q
            
            # Spawn a number of sampling trajectories using this equilibrated position as the starting point
            logging.info('Spawning {0:d} sampling trajectories at xi = {1:g}...'.format(numberOfTrajectories, xi_current))
            args = (self, xi_current, q, equilibrationSteps, evolutionSteps, saveTrajectories)
            for trajectory in range(numberOfTrajectories):
                results.append(pool.apply_async(runUmbrellaTrajectory, args))           

            logging.info('')
            
        # Now start at xi = 1 and move in the xi < 1 direction, using the
        # result of the previous xi as the initial position for the next xi
        q = q_initial[:,:,:,start]
        for l in range(start - 1, -1, -1):
            xi_current = xi_list[l]
            
            # Equilibrate in this window
            logging.info('Equilibrating trajectory at xi = {0:g} for {1:g} ps...'.format(xi_current, equilibrationSteps * self.dt * 2.418884326505e-5))
            p = self.sampleMomentum()
            result = system.equilibrate(0, p, q, equilibrationSteps, xi_current, self.potential, False, saveTrajectories)
            logging.info('Finished equilibrating trajectory at xi = {0:g}.'.format(xi_current))
            q_initial[:,:,:,l] = q

            # Spawn a number of sampling trajectories using this equilibrated position as the starting point
            logging.info('Spawning {0:d} sampling trajectories at xi = {1:g}...'.format(numberOfTrajectories, xi_current))
            args = (self, xi_current, q, equilibrationSteps, evolutionSteps, saveTrajectories)
            for trajectory in range(numberOfTrajectories):
                results.append(pool.apply_async(runUmbrellaTrajectory, args))           

            logging.info('')
        
        # Wait for each trajectory to finish, then update the mean and variance
        count = 0
        f = open('reaction_coordinate.dat', 'w')
        indices = range(start, Nxi)
        indices.extend(range(start - 1, -1, -1))
        for l in indices:
            xi_current = xi_list[l]
            logging.info('Processing {0:d} trajectories at xi = {1:g}...'.format(numberOfTrajectories, xi_current))
            for trajectory in range(numberOfTrajectories):
                
                # This line will block until the trajectory finishes
                dav, dav2 = results[count].get()
                
                # Update the mean and variance with the results from this trajectory
                # Note that these are counted at each time step in each trajectory
                av[l] += dav
                av2[l] += dav2
                
                # Print the updated mean and variance to the log file
                av_temp = av[l] / ((trajectory+1) * evolutionSteps)
                av2_temp = av2[l] / ((trajectory+1) * evolutionSteps)
                logging.info('{0:7d} {1:15.5e} {2:15.5e} {3:15.5e}'.format(trajectory+1, av_temp, av2_temp, av2_temp - av_temp * av_temp))
    
                count += 1
                
            logging.info('Finished processing trajectories at xi = {0:g}...'.format(xi_current))
            
            f.write('{0:9.5f} {1:15.5e} {2:15.5e}\n'.format(xi_current, av_temp, av2_temp - av_temp * av_temp))
                
        f.close()
        
    def computeTransmissionCoefficient(self, T, Nbeads, dt, 
                                       equilibrationTime,
                                       xi_current,
                                       parentEvolutionTime,
                                       childrenPerSampling,
                                       childEvolutionTime,
                                       childSamplingTime,
                                       processes=1,
                                       saveParentTrajectory=False, 
                                       saveChildTrajectories=False):
        """
        Return the transmission coefficient by the recrossing method. In this
        approach, a constrained RPMD simulation is initiated in the presence of
        an Andersen thermostat to generate a series of independent 
        configurations with centroids on the transition state dividing surface
        :math:`\\bar{s}_1(\\mathbf{q}) = 0`. For each of these "parent"
        configurations, a set of "child" trajectories is spawned by sampling
        from a Maxwell-Boltzmann distribution, with each trajectory evolving
        in time without the dividing surface constraint or Andersen thermostat.
        The transmission coefficient is then computed via
        
        .. math::
        
            \\kappa^{(n)}(s_1) = \\lim_{t \\rightarrow \\infty} 
                \\left< \\bar{f}_{s_1}(\\mathbf{q})^{-1} \\bar{v}_{s_1}(\\mathbf{p}, \\mathbf{q}) h \\left[ \\bar{s}_1(\\mathbf{q}_t) \\right] \\right>
        
        where the bracketed quantity is averaged over a large number of the
        child trajectories spawned from a large number of parent configurations.
        
        The `saveParentTrajectory` and `saveChildTrajectories` flags enable
        saving of the parent trajectory and/or a sampling of child trajectories
        as XYZ data files for later visualization in programs such as VMD.
        This is off by default because it is very slow. 
        """
        
        # Set the parameters for the RPMD calculation
        self.beta = 4.35974417e-18 / (constants.kB * T)
        self.dt = dt / 2.418884326505e-5
        self.Nbeads = Nbeads
        self.kforce = 0.0
        equilibrationTime /= 2.418884326505e-5
        parentEvolutionTime /= 2.418884326505e-5
        childEvolutionTime /= 2.418884326505e-5
        childSamplingTime /= 2.418884326505e-5
        geometry = self.transitionState.geometry[:,:,0]
        self.xi_current = xi_current
        self.mode = 2
        
        # Create a pool of subprocesses to farm out the individual trajectories to
        pool = multiprocessing.Pool(processes=processes)
        results = []

        logging.info('*****************************')
        logging.info('RPMD transmission coefficient')
        logging.info('*****************************')
        logging.info('')
        
        equilibrationSteps = int(round(equilibrationTime / self.dt))
        parentEvolutionSteps = int(round(parentEvolutionTime / self.dt))
        childEvolutionSteps = int(round(childEvolutionTime / self.dt))
        childSamplingSteps = int(round(childSamplingTime / self.dt))
        
        logging.info('Parameters')
        logging.info('==========')
        logging.info('Temperature                             = {0:g} K'.format(T))
        logging.info('Number of beads                         = {0:d}'.format(Nbeads))
        logging.info('Reaction coordinate                     = {0:g}'.format(xi_current))
        logging.info('Time step                               = {0:g} ps'.format(self.dt * 2.418884326505e-5))
        logging.info('Length of parent trajectory             = {0:g} ps ({1:d} steps)'.format(parentEvolutionSteps * self.dt * 2.418884326505e-5, parentEvolutionSteps))
        logging.info('Initial parent equilibration time       = {0:g} ps ({1:d} steps)'.format(equilibrationSteps * self.dt * 2.418884326505e-5, equilibrationSteps))
        logging.info('Frequency of child trajectory sampling  = {0:g} ps ({1:d} steps)'.format(childSamplingSteps * self.dt * 2.418884326505e-5, childSamplingSteps))
        logging.info('Length of child trajectories            = {0:g} ps ({1:d} steps)'.format(childEvolutionSteps * self.dt * 2.418884326505e-5, childEvolutionSteps))
        logging.info('Number of children per sampling         = {0:d}'.format(childrenPerSampling))
        logging.info('')
        
        self.activate()
        
        # Generate initial position using transition state geometry
        # (All beads start at same position)
        q = numpy.zeros((3,self.Natoms,self.Nbeads), order='F')
        for i in range(3):
            for j in range(self.Natoms):
                for k in range(self.Nbeads):
                    q[i,j,k] = geometry[i,j]
        # Sample initial momentum from normal distribution
        p = self.sampleMomentum()
        
        # Equilibrate parent trajectory while constraining to dividing surface
        # and sampling from Andersen thermostat
        logging.info('Equilibrating parent trajectory for {0:g} ps...'.format(equilibrationSteps * self.dt * 2.418884326505e-5))
        result = system.equilibrate(0, p, q, equilibrationSteps, self.xi_current, self.potential, True, saveParentTrajectory)
        
        logging.info('Finished equilibrating parent trajectory.')
        logging.info('')
        
        # Initialize parameters used to compute recrossing factor
        kappa_num = numpy.zeros(childEvolutionSteps, order='F')
        kappa_denom = numpy.array(0.0, order='F')
        
        # Continue evolving parent trajectory, interrupting to sample sets of
        # child trajectories in order to update the recrossing factor
        for iter in range(parentEvolutionSteps / childSamplingSteps + 1):
            
            logging.info('Sampling {0} child trajectories at {1:g} ps...'.format(childrenPerSampling, iter * childSamplingSteps * self.dt * 2.418884326505e-5))

            # Sample a number of child trajectories using the current parent
            # configuration
            results = []
            saveChildTrajectory = saveChildTrajectories
            for childCount in range(childrenPerSampling / 2):
                q_child = numpy.array(q.copy(), order='F')
                p_child = self.sampleMomentum()
                
                args = (self, xi_current, -p_child, q_child, childEvolutionSteps, saveChildTrajectory)
                results.append(pool.apply_async(runRecrossingTrajectory, args))           

                saveChildTrajectory = False
                
                args = (self, xi_current, p_child, q_child, childEvolutionSteps, saveChildTrajectory)
                results.append(pool.apply_async(runRecrossingTrajectory, args))           

            for childCount in range(childrenPerSampling):
                # This line will block until the child trajectory finishes
                num, denom = results[childCount].get()
                # Update the numerator and denominator of the recrossing factor expression
                kappa_num += num
                kappa_denom += denom
        
            logging.info('Finished sampling {0} child trajectories at {1:g} ps.'.format(childrenPerSampling, iter * childSamplingSteps * self.dt * 2.418884326505e-5))
            
            f = open('recrossing_factor.dat', 'w')
            for childStep in range(childEvolutionSteps):
                f.write('{0:11.3f} {1:11.6f} {2:11.6f}\n'.format(
                    childStep * self.dt * 2.418884326505e-2,
                    kappa_num[childStep] / kappa_denom,
                    kappa_num[childStep] / ((iter+1) * childrenPerSampling),
                ))
            f.close()
            
            logging.info('Current value of transmission coefficient = {0:.6f}'.format(kappa_num[-1] / kappa_denom))
            logging.info('')
            
            # Further evolve parent trajectory while constraining to dividing
            # surface and sampling from Andersen thermostat
            logging.info('Evolving parent trajectory to {0:g} ps...'.format((iter+1) * childSamplingSteps * self.dt * 2.418884326505e-5))
            result = system.equilibrate(0, p, q, childSamplingSteps, self.xi_current, self.potential, True, saveParentTrajectory)
        
        logging.info('Finished evolving parent trajectory for {0:g} ps...'.format(parentEvolutionSteps * self.dt * 2.418884326505e-5))
        logging.info('')
        
        logging.info('Result of recrossing factor calculation:')
        logging.info('')
        logging.info('=========== =========== ===========')
        logging.info('Time (fs)   kappa (new) kappa (old)')
        logging.info('=========== =========== ===========')
        
        for childStep in range(childEvolutionSteps):
            logging.info('{0:11.3f} {1:11.6f} {2:11.6f}'.format(
                childStep * self.dt * 2.418884326505e-2,
                kappa_num[childStep] / kappa_denom,
                kappa_num[childStep] / ((iter+1) * childrenPerSampling),
            ))
        logging.info('=========== =========== ===========')
        logging.info('')

        logging.info('Final value of transmission coefficient = {0:.6f}'.format(kappa_num[-1] / kappa_denom))
        logging.info('')
        
        return kappa_num[-1] / kappa_denom
    
    def sampleMomentum(self):
        """
        Return a pseudo-random sampling of momenta from a Boltzmann 
        distribution at the temperature of interest.
        """
        return system.sample_momentum(self.mass, self.beta, self.Nbeads)
