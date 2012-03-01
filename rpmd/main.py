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

import sys
import math
import numpy
import logging
import multiprocessing

import rpmd.constants as constants
import rpmd.quantity as quantity

from rpmd._main import *

################################################################################

class RPMDError(Exception):
    """
    An exception raised when an error occurs during an RPMD simulation. Pass a
    string describing the circumstances of the exceptional behavior.
    """
    pass

################################################################################

def runUmbrellaTrajectory(rpmd, xi_current, q, equilibrationSteps, evolutionSteps, kforce, saveTrajectory):
    """
    Run an individual umbrella integration trajectory, returning the sum of the
    first and second moments of the reaction coordinate at each time step.
    """
    rpmd.activate()
    p = rpmd.sampleMomentum()
    result = system.equilibrate(0, p, q, equilibrationSteps, xi_current, rpmd.potential, kforce, False, saveTrajectory)
    dav, dav2, result = system.umbrella_trajectory(0, p, q, evolutionSteps, xi_current, rpmd.potential, kforce, saveTrajectory)
    return dav, dav2, evolutionSteps

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

################################################################################

class Window:
    """
    A representation of a window along the reaction coordinate used for
    umbrella sampling. The attributes are:
    
    =========================== ================================================
    Attribute                   Description
    =========================== ================================================
    `xi`                        The value of the reaction coordinate in the center of the window
    `kforce`                    The umbrella integration force constant for this window
    `trajectories`              The number of independent sampling trajectories to run for this window
    `equilibrationTime`         The equilibration time (no sampling) in each trajectory
    `evolutionTime`             The evolution time (with sampling) in each trajectory
    --------------------------- ------------------------------------------------
    `count`                     The number of samples taken
    `av`                        The mean of the reaction coordinate times the number of samples
    `av2`                       The variance of the reaction coordinate times the number of samples
    =========================== ================================================    
    
    """
    
    def __init__(self, xi, kforce, trajectories, equilibrationTime, evolutionTime):
        # These parameters control the umbrella sampling trajectories
        self.xi = xi
        self.kforce = kforce
        self.trajectories = trajectories
        self.equilibrationTime = float(quantity.convertTime(equilibrationTime, "ps")) / 2.418884326505e-5
        self.evolutionTime = float(quantity.convertTime(evolutionTime, "ps")) / 2.418884326505e-5
        # The parameters store the results of the sampling
        self.count = 0
        self.av = 0.0
        self.av2 = 0.0

################################################################################

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
    `transitionStates`          The dividing surface(s) near the transition state, as a list of :class:`TransitionState` objects
    `potential`                 A function that computes the potential and forces for a given position
    --------------------------- ------------------------------------------------
    `beta`                      The reciprocal temperature of the RPMD simulation
    `dt`                        The time step to use in the RPMD simulation
    `Nbeads`                    The number of beads per atom in the RPMD simulation
    `xi_current`                The current value of the reaction coordinate
    `mode`                      A flag indicating the type of RPMD calculation currently underway (1 = umbrella, 2 = recrossing)
    =========================== ================================================
    
    """

    def __init__(self, label, T, reactants, transitionStates, potential):
        """
        Initialize an RPMD object. The `mass` of each atom should be given in
        g/mol, while the `Rinf` value should be given in angstroms. (They will
        be converted to atomic units.)
        """
        self.label = label
        self.T = float(quantity.convertTemperature(T, "K"))
        self.mass = reactants.mass
        self.Natoms = len(self.mass)
        self.reactants = reactants
        self.transitionStates = transitionStates or []
        self.potential = potential
        
        self.beta = 4.35974417e-18 / (constants.kB * self.T)
        self.dt = 0
        self.Nbeads = 0
        self.xi_current = 0
        self.mode = 0
        
        self.umbrellaConfigurations = None
        self.umbrellaWindows = None
        self.potentialOfMeanForce = None
        self.recrossingFactor = None
        
        self.initializeLog()
        
    def initializeLog(self, verbose=logging.INFO):
        """
        Set up a logger for RPMD to use to print output to stdout. The
        `verbose` parameter is an integer specifying the amount of log text seen
        at the console; the levels correspond to those of the :data:`logging` module.
        """
        # Create logger
        logger = logging.getLogger()
        logger.setLevel(verbose)
    
        # Create console handler; send everything to stdout rather than stderr
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(verbose)
    
        logging.addLevelName(logging.CRITICAL, 'Critical: ')
        logging.addLevelName(logging.ERROR, 'Error: ')
        logging.addLevelName(logging.WARNING, 'Warning: ')
        logging.addLevelName(logging.INFO, '')
        logging.addLevelName(logging.DEBUG, '')
        logging.addLevelName(0, '')
    
        # Create formatter and add to handlers
        formatter = logging.Formatter('%(levelname)s%(message)s')
        ch.setFormatter(formatter)
        
        # Remove old handlers before adding ours
        while logger.handlers:
            logger.removeHandler(logger.handlers[0])
    
        # Add handlers to logger
        logger.addHandler(ch)

    def activate(self):
        """
        Set this object as the active RPMD system in the Fortran layer. Note
        that the dividing surface properties must be set in the corresponding
        modules in ``rpmd._main``, *not* in ``rpmd._surface``.
        """
        Natoms = self.mass.shape[0]
        system.dt = self.dt
        system.beta = self.beta
        system.mass[0:Natoms] = self.mass
        system.mode = self.mode
        self.reactants.activate(module=reactants)
        self.thermostat.activate(module=system)
        
        Nts = len(self.transitionStates)
        Nforming_bonds = max([ts.formingBonds.shape[1] for ts in self.transitionStates])
        Nbreaking_bonds = max([ts.breakingBonds.shape[1] for ts in self.transitionStates])

        formingBonds = numpy.zeros((Nts,Nforming_bonds,2))
        breakingBonds = numpy.zeros((Nts,Nbreaking_bonds,2))
        formingBondLengths = numpy.zeros((Nts,Nforming_bonds))
        breakingBondLengths = numpy.zeros((Nts,Nbreaking_bonds))
        
        for n, ts in enumerate(self.transitionStates):
            formingBonds[n,:,:] = ts.formingBonds
            breakingBonds[n,:,:] = ts.breakingBonds
            formingBondLengths[n,:] = ts.formingBondLengths
            breakingBondLengths[n,:] = ts.breakingBondLengths
        
        transition_state.number_of_transition_states = Nts
        transition_state.number_of_forming_bonds = Nforming_bonds
        transition_state.forming_bonds[0:Nts,0:Nforming_bonds,:] = formingBonds
        transition_state.forming_bond_lengths[0:Nts,0:Nforming_bonds] = formingBondLengths
        transition_state.number_of_breaking_bonds = Nbreaking_bonds
        transition_state.breaking_bonds[0:Nts,0:Nbreaking_bonds,:] = breakingBonds
        transition_state.breaking_bond_lengths[0:Nts,0:Nbreaking_bonds] = breakingBondLengths
    
    def generateUmbrellaConfigurations(self, 
                                       dt, 
                                       evolutionTime,
                                       xi_list,
                                       kforce,
                                       thermostat):
        """
        Generate a set of configurations along the reaction coordinate for
        future use in RPMD umbrella sampling. The algorithm starts near the
        transition state dividing surface and moves away in either direction,
        running an RPMD equilibration in each window to obtain the appropriate
        configuration in that window. That configuration is then used as the
        initial position for determining the configuration in the next window.
        """
        
        dt = float(quantity.convertTime(dt, "ps")) / 2.418884326505e-5
        evolutionTime = float(quantity.convertTime(evolutionTime, "ps")) / 2.418884326505e-5
        
        # Set the parameters for the RPMD calculation
        self.dt = dt
        self.Nbeads = 1
        xi_list = numpy.array(xi_list)
        Nxi = len(xi_list)
        geometry = self.transitionStates[0].geometry
        self.thermostat = thermostat
        self.mode = 1
        
        if isinstance(kforce, float):
            kforce = numpy.ones_like(xi_list) * kforce
        
        logging.info('****************************')
        logging.info('RPMD umbrella configurations')
        logging.info('****************************')
        logging.info('')
        
        evolutionSteps = int(round(evolutionTime / self.dt))
        
        logging.info('Parameters')
        logging.info('==========')
        logging.info('Temperature                             = {0:g} K'.format(self.T))
        logging.info('Number of beads                         = {0:d}'.format(self.Nbeads))
        logging.info('Time step                               = {0:g} ps'.format(self.dt * 2.418884326505e-5))
        logging.info('Number of umbrella integration windows  = {0:d}'.format(Nxi))
        logging.info('Trajectory evolution time               = {0:g} ps ({1:d} steps)'.format(evolutionSteps * self.dt * 2.418884326505e-5, evolutionSteps))
        logging.info('')

        # Only use one bead to generate initial positions in each window
        # (We will equilibrate within each window to allow the beads to separate)
        self.Nbeads = 1
        self.activate()

        # Generate initial position using transition state geometry
        # (All beads start at same position)
        q = numpy.zeros((3,self.Natoms,self.Nbeads), order='F')
        for i in range(3):
            for j in range(self.Natoms):
                for k in range(self.Nbeads):
                    q[i,j,k] = geometry[i,j]

        # Find the window nearest to the transition state dividing surface
        for start in range(Nxi):
            if xi_list[start] >= 1:
                break
        
        # Equilibrate in each window to determine the initial positions
        # First start at xi = 1 and move in the xi > 1 direction, using the
        # result of the previous xi as the initial position for the next xi
        q_initial = numpy.zeros((3,self.Natoms,Nxi), order='F')
        for l in range(start, Nxi):
            xi_current = xi_list[l]
            
            # Equilibrate in this window
            logging.info('Generating configuration at xi = {0:g} for {1:g} ps...'.format(xi_current, evolutionSteps * self.dt * 2.418884326505e-5))
            p = self.sampleMomentum()
            result = system.equilibrate(0, p, q, evolutionSteps, xi_current, self.potential, kforce[l], False, False)
            logging.info('Finished generating configuration at xi = {0:g}.'.format(xi_current))
            q_initial[:,:,l] = q[:,:,0]
                        
        # Now start at xi = 1 and move in the xi < 1 direction, using the
        # result of the previous xi as the initial position for the next xi
        q[:,:,0] = q_initial[:,:,start]
        for l in range(start - 1, -1, -1):
            xi_current = xi_list[l]
            
            # Equilibrate in this window
            logging.info('Generating configuration at xi = {0:g} for {1:g} ps...'.format(xi_current, evolutionSteps * self.dt * 2.418884326505e-5))
            p = self.sampleMomentum()
            result = system.equilibrate(0, p, q, evolutionSteps, xi_current, self.potential, kforce[l], False, False)
            logging.info('Finished generating configuration at xi = {0:g}.'.format(xi_current))
            q_initial[:,:,l] = q[:,:,0]
        
        # Store the computed configurations on the object for future use in
        # umbrella sampling
        self.umbrellaConfigurations = []
        for l in range(Nxi):
            xi_current = xi_list[l]
            q_current = self.cleanGeometry(q_initial[:,:,l])
            self.umbrellaConfigurations.append((xi_current, q_current))
            logging.info('Configuration at xi = {0:g}:'.format(xi_current))
            for j in range(self.Natoms):
                logging.info('{0:5} {1:11.6f} {2:11.6f} {3:11.6f}'.format(self.reactants.atoms[j], q_current[0,j], q_current[1,j], q_current[2,j]))                
            logging.info('')
        
    def conductUmbrellaSampling(self, 
                                Nbeads, 
                                dt, 
                                windows,
                                thermostat,
                                processes=1,
                                saveTrajectories=False):
        """
        Return the value of the static factor :math:`p^{(n)}(s_1, s_0)` as
        computed using umbrella integration.
        """
        
        # Don't continue if the user hasn't generated the initial configurations yet
        if not self.umbrellaConfigurations:
            raise RPMDError('You must run generateUmbrellaConfigurations() before running computeStaticFactor().')
        
        # Set the parameters for the RPMD calculation
        self.dt = dt = float(quantity.convertTime(dt, "ps")) / 2.418884326505e-5
        self.Nbeads = Nbeads
        Nwindows = len(windows)
        self.thermostat = thermostat
        self.mode = 1
        
        # Create a pool of subprocesses to farm out the individual trajectories to
        pool = multiprocessing.Pool(processes=processes)
        results = []

        logging.info('******************')
        logging.info('RPMD static factor')
        logging.info('******************')
        logging.info('')
        
        logging.info('Parameters')
        logging.info('==========')
        logging.info('Temperature                             = {0:g} K'.format(self.T))
        logging.info('Number of beads                         = {0:d}'.format(Nbeads))
        logging.info('Time step                               = {0:g} ps'.format(self.dt * 2.418884326505e-5))
        logging.info('Number of umbrella integration windows  = {0:d}'.format(Nwindows))
        logging.info('')

        self.activate()

        self.umbrellaWindows = windows

        # Spawn a number of sampling trajectories in each window
        for window in windows:
            
            equilibrationSteps = int(round(window.equilibrationTime / self.dt))
            evolutionSteps = int(round(window.evolutionTime / self.dt))
        
            # Load initial configuration using results from generateUmbrellaConfigurations()
            q = numpy.empty((3,self.Natoms,self.Nbeads), order='F')
            for xi, q_initial in self.umbrellaConfigurations:
                if xi >= window.xi:
                    break
            for k in range(self.Nbeads):
                q[:,:,k] = q_initial
            
            logging.info('Spawning {0:d} sampling trajectories at xi = {1:g}...'.format(window.trajectories, window.xi))
            args = (self, window.xi, q, equilibrationSteps, evolutionSteps, window.kforce, saveTrajectories)
            for trajectory in range(window.trajectories):
                results.append(pool.apply_async(runUmbrellaTrajectory, args))           

        logging.info('')
                    
        # Wait for each trajectory to finish, then update the mean and variance
        count = 0
        f = open('reaction_coordinate.dat', 'w')
        for window in windows:
            logging.info('Processing {0:d} trajectories at xi = {1:g}...'.format(window.trajectories, window.xi))
            for trajectory in range(window.trajectories):
                
                # This line will block until the trajectory finishes
                dav, dav2, dcount = results[count].get()
                
                # Update the mean and variance with the results from this trajectory
                # Note that these are counted at each time step in each trajectory
                window.av += dav
                window.av2 += dav2
                window.count += dcount
                
                # Print the updated mean and variance to the log file
                av_temp = window.av / window.count
                av2_temp = window.av2 / window.count
                logging.info('{0:5d} {1:11d} {2:15.8f} {3:15.8f} {4:15.5e}'.format(trajectory+1, window.count, av_temp, av2_temp, av2_temp - av_temp * av_temp))
    
                count += 1
                
            logging.info('Finished processing trajectories at xi = {0:g}...'.format(window.xi))
            
            f.write('{0:9.5f} {1:15.5e} {2:15.5e}\n'.format(window.xi, av_temp, av2_temp - av_temp * av_temp))
                
        f.close()
        
        logging.info('')
        
    def computePotentialOfMeanForce(self, xi_min, xi_max, bins):
        """
        Compute the potential of mean force of the system at the given
        temperature by integrating over the given reaction coordinate range
        using the given number of bins. This requires that you have previously
        used umbrella sampling to determine the mean and variance in each bin.
        """
        # Don't continue if the user hasn't generated the umbrella sampling yet
        if not self.umbrellaWindows:
            raise RPMDError('You must run conductUmbrellaSampling() before running computePotentialOfMeanForce().')
        
        logging.info('****************************')
        logging.info('RPMD potential of mean force')
        logging.info('****************************')
        logging.info('')
        
        logging.info('Parameters')
        logging.info('==========')
        logging.info('Temperature                             = {0:g} K'.format(self.T))
        logging.info('Lower bound of reaction coordinate      = {0:g}'.format(xi_min))
        logging.info('Upper bound of reaction coordinate      = {0:g}'.format(xi_max))
        logging.info('Number of bins                          = {0:d}'.format(bins))
        logging.info('')
        
        Nwindows = len(self.umbrellaWindows)
        
        xi_list = numpy.linspace(xi_min, xi_max, bins, True)
        
        # Count the number of bins in each window
        N = numpy.zeros(Nwindows)
        for l in range(1, Nwindows-1):
            xi_left = 0.5 * (self.umbrellaWindows[l-1].xi + self.umbrellaWindows[l].xi)
            xi_right = 0.5 * (self.umbrellaWindows[l].xi + self.umbrellaWindows[l+1].xi)
            N[l] = sum([1 for xi in xi_list if xi_left <= xi < xi_right])
        xi_right = 0.5 * (self.umbrellaWindows[0].xi + self.umbrellaWindows[1].xi)
        N[0] = sum([1 for xi in xi_list if xi < xi_right])
        xi_left = 0.5 * (self.umbrellaWindows[-1].xi + self.umbrellaWindows[-2].xi)
        N[-1] = sum([1 for xi in xi_list if xi_left <= xi])
        
        # Compute the slope in each bin
        dA = numpy.zeros(bins)
        p = numpy.zeros(Nwindows)
        dA0 = numpy.zeros(Nwindows)
        for n, xi in enumerate(xi_list):
            for l, window in enumerate(self.umbrellaWindows):
                xi_mean = window.av / window.count
                xi_var = window.av2 / window.count
                xi_window = window.xi
                kforce = window.kforce
                p[l] = 1.0 / numpy.sqrt(2 * constants.pi * xi_var) * numpy.exp(-0.5 * (xi - xi_mean)**2 / xi_var) 
                dA0[l] = (1.0 / self.beta) * (xi - xi_mean) / xi_var - kforce * (xi - xi_window)
            dA[n] = numpy.sum(N * p * dA0) / numpy.sum(N * p)
            
        # Now integrate numerically to get the potential of mean force
        self.potentialOfMeanForce = numpy.zeros((2,bins))
        for n, xi in enumerate(xi_list):
            self.potentialOfMeanForce[0,n] = xi_list[n]
            self.potentialOfMeanForce[1,n] = numpy.trapz(dA[:n], xi_list[:n])
             
        logging.info('Result of potential of mean force calculation:')
        logging.info('')
        logging.info('=========== ===========')
        logging.info('Rxn coord   PMF (eV)')
        logging.info('=========== ===========')
        for n in range(0, self.potentialOfMeanForce.shape[1], 10):
            logging.info('{0:11.6f} {1:11.6f}'.format(
                self.potentialOfMeanForce[0,n],
                self.potentialOfMeanForce[1,n] * 27.211,
            ))
        logging.info('=========== ===========')
        logging.info('')

    def computeTransmissionCoefficient(self, 
                                       Nbeads, 
                                       dt, 
                                       equilibrationTime,
                                       childTrajectories,
                                       childrenPerSampling,
                                       childEvolutionTime,
                                       childSamplingTime,
                                       thermostat,
                                       processes=1,
                                       xi_current=None,
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
        
        # If xi_current not specified, use the maximum of the potential of mean force
        if xi_current is None:
            if self.potentialOfMeanForce is None:
                raise RPMDError('You must run computePotentialOfMeanForce() or specify xi_current before running computeTransmissionCoefficient().')
            index = numpy.argmax(self.potentialOfMeanForce[1,:])
            xi_current = self.potentialOfMeanForce[0,index]
        
        dt = float(quantity.convertTime(dt, "ps")) / 2.418884326505e-5
        equilibrationTime = float(quantity.convertTime(equilibrationTime, "ps")) / 2.418884326505e-5
        childEvolutionTime = float(quantity.convertTime(childEvolutionTime, "ps")) / 2.418884326505e-5
        childSamplingTime = float(quantity.convertTime(childSamplingTime, "ps")) / 2.418884326505e-5

        # Set the parameters for the RPMD calculation
        self.dt = dt
        self.Nbeads = Nbeads
        self.kforce = 0.0
        geometry = self.transitionStates[0].geometry
        self.xi_current = xi_current
        self.thermostat = thermostat
        self.mode = 2
        
        # Create a pool of subprocesses to farm out the individual trajectories to
        pool = multiprocessing.Pool(processes=processes)
        results = []

        logging.info('*****************************')
        logging.info('RPMD transmission coefficient')
        logging.info('*****************************')
        logging.info('')
        
        equilibrationSteps = int(round(equilibrationTime / self.dt))
        childEvolutionSteps = int(round(childEvolutionTime / self.dt))
        childSamplingSteps = int(round(childSamplingTime / self.dt))
        
        logging.info('Parameters')
        logging.info('==========')
        logging.info('Temperature                             = {0:g} K'.format(self.T))
        logging.info('Number of beads                         = {0:d}'.format(Nbeads))
        logging.info('Reaction coordinate                     = {0:g}'.format(xi_current))
        logging.info('Time step                               = {0:g} ps'.format(self.dt * 2.418884326505e-5))
        logging.info('Total number of child trajectories      = {0:d}'.format(childTrajectories))
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
        result = system.equilibrate(0, p, q, equilibrationSteps, self.xi_current, self.potential, 0.0, True, saveParentTrajectory)
        
        logging.info('Finished equilibrating parent trajectory.')
        logging.info('')
        
        # Initialize parameters used to compute recrossing factor
        kappa_num = numpy.zeros(childEvolutionSteps, order='F')
        kappa_denom = numpy.array(0.0, order='F')
        
        # Continue evolving parent trajectory, interrupting to sample sets of
        # child trajectories in order to update the recrossing factor
        childCount = 0; parentIter = 0
        while childCount < childTrajectories:
            
            logging.info('Sampling {0} child trajectories at {1:g} ps...'.format(childrenPerSampling, parentIter * childSamplingSteps * self.dt * 2.418884326505e-5))

            # Sample a number of child trajectories using the current parent
            # configuration
            results = []
            saveChildTrajectory = saveChildTrajectories
            for child in range(childrenPerSampling / 2):
                q_child = numpy.array(q.copy(), order='F')
                p_child = self.sampleMomentum()
                
                args = (self, xi_current, -p_child, q_child, childEvolutionSteps, saveChildTrajectory)
                results.append(pool.apply_async(runRecrossingTrajectory, args))           
                childCount += 1
                
                saveChildTrajectory = False
                
                args = (self, xi_current, p_child, q_child, childEvolutionSteps, saveChildTrajectory)
                results.append(pool.apply_async(runRecrossingTrajectory, args))
                childCount += 1         

            for child in range(childrenPerSampling):
                # This line will block until the child trajectory finishes
                num, denom = results[child].get()
                # Update the numerator and denominator of the recrossing factor expression
                kappa_num += num
                kappa_denom += denom
        
            logging.info('Finished sampling {0} child trajectories at {1:g} ps.'.format(childrenPerSampling, parentIter * childSamplingSteps * self.dt * 2.418884326505e-5))
            
            f = open('recrossing_factor.dat', 'w')
            for childStep in range(childEvolutionSteps):
                f.write('{0:11.3f} {1:11.6f} {2:11.6f}\n'.format(
                    childStep * self.dt * 2.418884326505e-2,
                    kappa_num[childStep] / kappa_denom,
                    kappa_num[childStep] / ((parentIter+1) * childrenPerSampling),
                ))
            f.close()
            
            logging.info('Current value of transmission coefficient = {0:.6f}'.format(kappa_num[-1] / kappa_denom))
            logging.info('')
            
            # Further evolve parent trajectory while constraining to dividing
            # surface and sampling from Andersen thermostat
            logging.info('Evolving parent trajectory to {0:g} ps...'.format((parentIter+1) * childSamplingSteps * self.dt * 2.418884326505e-5))
            result = system.equilibrate(0, p, q, childSamplingSteps, self.xi_current, self.potential, 0.0, True, saveParentTrajectory)
        
            parentIter += 1
        
        logging.info('Finished sampling of {0:d} child trajectories.'.format(childTrajectories))
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
                kappa_num[childStep] / ((parentIter+1) * childrenPerSampling),
            ))
        logging.info('=========== =========== ===========')
        logging.info('')

        logging.info('Final value of transmission coefficient = {0:.6f}'.format(kappa_num[-1] / kappa_denom))
        logging.info('')
        
        self.recrossingFactor = kappa_num[-1] / kappa_denom
        
        return self.recrossingFactor
    
    def computeRPMDRateCoefficient(self):
        """
        Compute the value of the RPMD rate coefficient.
        """
        
        if self.potentialOfMeanForce is None or self.recrossingFactor is None:
            raise RPMDError('You must run computePotentialOfMeanForce() and computeTransmissionCoefficient() before running computeRPMDRateCoefficient().')
        
        logging.info('*********************')
        logging.info('RPMD rate coefficient')
        logging.info('*********************')
        logging.info('')

        logging.info('Parameters')
        logging.info('==========')
        logging.info('Temperature                             = {0:g} K'.format(self.T))
        logging.info('')
        
        fromAtomicUnits = 1e6 / ((5.2917721092e-11)**3 / 2.418884326505e-17)
        
        # Compute the rate coefficient using the reactant dividing surface
        Rinf = self.reactants.Rinf
        mA = self.reactants.totalMass1
        mB = self.reactants.totalMass2
        mu = mA * mB / (mA + mB)
        k_QTST = 4 * constants.pi * Rinf * Rinf / numpy.sqrt(2 * constants.pi * self.beta * mu)
        
        logging.info('Result of RPMD rate coefficient calculation:')
        logging.info('')
        logging.info('k(T;s0) (QTST)                          = {0:g} cm^3/(mol*s)'.format(k_QTST * fromAtomicUnits))
        
        # Compute the static factor
        W1 = numpy.max(self.potentialOfMeanForce[1,:])
        W0 = self.potentialOfMeanForce[1,0]
        staticFactor = numpy.exp(-self.beta * (W1 - W0))
        
        # Correct the rate coefficient to the transition state dividing surface
        k_QTST *= staticFactor
        
        # Correct the rate coefficient for recrossings
        k_RPMD = k_QTST * self.recrossingFactor
        
        logging.info('Static factor                           = {0:g}'.format(staticFactor))
        logging.info('k(T;s1) (QTST)                          = {0:g} cm^3/(mol*s)'.format(k_QTST * fromAtomicUnits))
        logging.info('Dynamic (recrossing) factor             = {0:g}'.format(self.recrossingFactor))
        logging.info('k(T) (RPMD)                             = {0:g} cm^3/(mol*s)'.format(k_RPMD * fromAtomicUnits))
        logging.info('')
        
        logging.info('Final value of rate coefficient = {0:g} cm^3/(mol*s)'.format(k_RPMD * fromAtomicUnits))
        logging.info('')

        return k_RPMD
    
    def sampleMomentum(self):
        """
        Return a pseudo-random sampling of momenta from a Boltzmann 
        distribution at the temperature of interest.
        """
        return system.sample_momentum(self.mass, self.beta, self.Nbeads)

    def getCenterOfMass(self, position):
        """
        Return the center of mass for the given `position`.
        """
        cm = numpy.zeros(position.shape[0])
        mass = self.reactants.mass
        
        for i in range(position.shape[0]):
            for j in range(position.shape[1]):
                cm[i] += position[i,j] * mass[j]
                    
        cm /= numpy.sum(mass)
        
        return cm

    def cleanGeometry(self, geometry):
        """
        Return a copy of the geometry translated so that the center of mass
        is at the origin.
        """
        newGeometry = geometry.copy()
        cm = self.getCenterOfMass(geometry)
        for j in range(self.Natoms):
            newGeometry[:,j] -= cm
        return newGeometry
