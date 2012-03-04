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

import os
import os.path
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

def runUmbrellaTrajectory(rpmd, xi_current, p, q, equilibrationSteps, evolutionSteps, kforce, saveTrajectory):
    """
    Run an individual umbrella integration trajectory, returning the sum of the
    first and second moments of the reaction coordinate at each time step.
    """
    rpmd.activate()
    p = numpy.asfortranarray(p)
    q = numpy.asfortranarray(q)
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
    p = numpy.asfortranarray(p)
    q = numpy.asfortranarray(q)
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

    def __init__(self, label, T, reactants, transitionStates, potential, outputDirectory='.'):
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
        self.outputDirectory = os.path.abspath(outputDirectory)
        
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
        Nforming_bonds = max([ts.formingBonds.shape[0] for ts in self.transitionStates])
        Nbreaking_bonds = max([ts.breakingBonds.shape[0] for ts in self.transitionStates])

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

        # Set up output files and directory
        workingDirectory = self.createWorkingDirectory()
        configurationsFilename = os.path.join(workingDirectory, 'umbrella_configurations.dat')

        # Look for existing output file for this calculation
        # If a file exists, we won't repeat the calculation
        if os.path.exists(configurationsFilename):
            logging.info('Loading saved output from {0}'.format(configurationsFilename))
            xi_list0, q_initial0, equilibrationSteps0 = self.loadUmbrellaConfigurations(configurationsFilename)
            if xi_list.shape[0] == xi_list0.shape[0] and all(numpy.abs(xi_list - xi_list0) < 1e-6):
                logging.info('Using results of previously saved umbrella configurations.')
                logging.info('')
                xi_list = xi_list0
                q_initial = q_initial0
                self.umbrellaConfigurations = []
                for l in range(Nxi):
                    xi_current = xi_list[l]
                    q_current = self.cleanGeometry(q_initial[:,:,l])
                    self.umbrellaConfigurations.append((xi_current, q_current))
                return
            else:
                logging.info('NOT using results of previously saved umbrella configurations.')           
        else:
            logging.info('Output will be saved to {0}'.format(configurationsFilename))
        logging.info('')

        # Only use one bead to generate initial positions in each window
        # (We will equilibrate within each window to allow the beads to separate)
        self.Nbeads = 1
        self.activate()

        # Seed the random number generator
        random_init()

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
        
        self.saveUmbrellaConfigurations(configurationsFilename, evolutionSteps)
    
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

        logging.info('**********************')
        logging.info('RPMD umbrella sampling')
        logging.info('**********************')
        logging.info('')
        
        logging.info('Parameters')
        logging.info('==========')
        logging.info('Temperature                             = {0:g} K'.format(self.T))
        logging.info('Number of beads                         = {0:d}'.format(Nbeads))
        logging.info('Time step                               = {0:g} ps'.format(self.dt * 2.418884326505e-5))
        logging.info('Number of umbrella integration windows  = {0:d}'.format(Nwindows))
        logging.info('')

        # Set up output files and directory
        workingDirectory = self.createWorkingDirectory()
        umbrellaFilename = os.path.join(workingDirectory, 'umbrella_sampling_{0:d}.dat'.format(self.Nbeads))

        # Look for existing output file for this calculation
        # If a file exists, we won't repeat the calculation
        if os.path.exists(umbrellaFilename):
            logging.info('Loading saved output from {0}'.format(umbrellaFilename))
            xi_list0, av_list0, av2_list0, count_list0 = self.loadUmbrellaSampling(umbrellaFilename)
            Nxi0 = xi_list0.shape[0]
            for window in windows:
                for l in range(Nxi0):
                    if abs(xi_list0[l] - window.xi) < 1e-6:
                        window.av += av_list0[l]
                        window.av2 += av2_list0[l]
                        window.count += count_list0[l]
                        break
        else:
            logging.info('Output will be saved to {0}'.format(umbrellaFilename))
        logging.info('')
        
        self.activate()

        # Seed the random number generator
        random_init()

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
            
            # Spawn sampling trajectories in this window
            # In order to get better statistics (and therefore faster
            # convergence), we want to give each trajectory a unique initial
            # position and momentum
            # To do this, we spawn one trajectory after each equilibration
            # period, giving up a (small) bit of parallelization in the name of
            # better statistics in the sampling
            windowEvolutionSteps = evolutionSteps - int(numpy.ceil(float(window.count) / window.trajectories))
            windowEquilibrationSteps = equilibrationSteps
            if windowEvolutionSteps <= 0:
                logging.info('Already sampled enough trajectories at xi = {1:g}.'.format(window.trajectories, window.xi))
                windowEquilibrationSteps = 0
            else:
                logging.info('Spawning {0:d} sampling trajectories at xi = {1:g}...'.format(window.trajectories, window.xi))
            for trajectory in range(window.trajectories):
                p = self.sampleMomentum()
                result = system.equilibrate(0, p, q, windowEquilibrationSteps, window.xi, self.potential, window.kforce, False, False)
                args = (self, window.xi, p, q, 0, windowEvolutionSteps, window.kforce, saveTrajectories)
                results.append(pool.apply_async(runUmbrellaTrajectory, args))           

        logging.info('')
                 
        # Wait for each trajectory to finish, then update the mean and variance
        count = 0
        f = open(umbrellaFilename, 'w')

        f.write('**********************\n')
        f.write('RPMD umbrella sampling\n')
        f.write('**********************\n\n')

        f.write('Temperature                             = {0:g} K\n'.format(self.T))
        f.write('Number of beads                         = {0:d}\n'.format(self.Nbeads))
        f.write('Time step                               = {0:g} ps\n'.format(self.dt * 2.418884326505e-5))
        f.write('Number of umbrella integration windows  = {0:d}\n\n'.format(len(self.umbrellaConfigurations)))
        
        f.write('========= =============== =============== =========== ============= =============\n')
        f.write('xi        total av        total av2       count       xi_mean       xi_var\n')
        f.write('========= =============== =============== =========== ============= =============\n')

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
                if dcount > 0:
                    logging.info('{0:5d} {1:11d} {2:15.8f} {3:15.8f} {4:15.5e}'.format(trajectory+1, window.count, av_temp, av2_temp, av2_temp - av_temp * av_temp))
    
                count += 1
                
            logging.info('Finished processing trajectories at xi = {0:g}...'.format(window.xi))
            
            f.write('{0:9.5f} {1:15.8e} {2:15.8e} {3:11d} {4:13.5e} {5:13.5e}\n'.format(window.xi, window.av, window.av2, window.count, av_temp, av2_temp - av_temp * av_temp))

        f.write('========= =============== =============== =========== ============= =============\n')
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
        
        # Count the number of sampling points in each window
        N = numpy.zeros(Nwindows)
        for l in range(Nwindows):
            N[l] = self.umbrellaWindows[l].count
        
        # Compute the slope in each bin
        dA = numpy.zeros(bins)
        p = numpy.zeros(Nwindows)
        dA0 = numpy.zeros(Nwindows)
        for n, xi in enumerate(xi_list):
            for l, window in enumerate(self.umbrellaWindows):
                av = window.av / window.count
                av2 = window.av2 / window.count
                xi_mean = av
                xi_var = av2 - av * av
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

    def computeRecrossingFactor(self, 
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
        Return the recrossing factor for the RPMD system. A constrained RPMD
        simulation is initiated in the presence of a thermostat to generate a
        series of independent configurations with centroids on the transition
        state dividing surface :math:`\\bar{s}_1(\\mathbf{q}) = 0`. For each of
        these "parent" configurations, a set of "child" trajectories is spawned
        by sampling from a Maxwell-Boltzmann distribution, with each trajectory
        evolving in time without the dividing surface constraint or thermostat.
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
                raise RPMDError('You must run computePotentialOfMeanForce() or specify xi_current before running computeRecrossingFactor().')
            index = numpy.argmax(self.potentialOfMeanForce[1,:])
            xi_current = self.potentialOfMeanForce[0,index]
        
        dt = float(quantity.convertTime(dt, "ps")) / 2.418884326505e-5
        equilibrationTime = float(quantity.convertTime(equilibrationTime, "ps")) / 2.418884326505e-5
        childEvolutionTime = float(quantity.convertTime(childEvolutionTime, "ps")) / 2.418884326505e-5
        childSamplingTime = float(quantity.convertTime(childSamplingTime, "ps")) / 2.418884326505e-5
        equilibrationSteps = int(round(equilibrationTime / dt))
        childEvolutionSteps = int(round(childEvolutionTime / dt))
        childSamplingSteps = int(round(childSamplingTime / dt))
        
        # Set the parameters for the RPMD calculation
        self.dt = dt
        self.Nbeads = Nbeads
        self.kforce = 0.0
        geometry = self.transitionStates[0].geometry
        self.xi_current = xi_current
        self.thermostat = thermostat
        self.mode = 2
        
        # Initialize parameters used to compute recrossing factor
        kappa_num = numpy.zeros(childEvolutionSteps, order='F')
        kappa_denom = numpy.array(0.0, order='F')
        childCount = 0
        
        # Create a pool of subprocesses to farm out the individual trajectories to
        pool = multiprocessing.Pool(processes=processes)
        results = []

        logging.info('**********************')
        logging.info('RPMD recrossing factor')
        logging.info('**********************')
        logging.info('')
        
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
        
        # Set up output files and directory
        workingDirectory = self.createWorkingDirectory()
        recrossingFilename = os.path.join(workingDirectory, 'recrossing_factor_{0:d}.dat'.format(Nbeads))

        # Look for existing output file for this calculation
        # If a file exists, we won't repeat the calculation unless more
        # child trajectories are requested
        if os.path.exists(recrossingFilename):
            logging.info('Loading saved output from {0}'.format(recrossingFilename))
            (T0, Nbeads0, xi_current0, dt0, kappa_num0, kappa_denom0, trajectoryCount0,
                childTrajectories0, equilibrationSteps0, childSamplingSteps0, 
                childEvolutionSteps0, childrenPerSampling0) = self.loadRecrossingFactor(recrossingFilename)
            if T0 == self.T and Nbeads0 == self.Nbeads and dt0 == self.dt and abs(xi_current0 - self.xi_current) < 1e-6:
                # We can use the old data
                logging.info('Including previously saved output in calculation.')
                kappa_num = kappa_num0
                kappa_denom = kappa_denom0
                childCount = trajectoryCount0
                logging.info('Saved output contained {0:d} child trajectories; {1:d} additional additional trajectories will be run.'.format(childCount, max(childTrajectories - childCount, 0)))           
            else:
                logging.info('NOT including previously saved output in calculation.')           
        else:
            logging.info('Output will be saved to {0}'.format(recrossingFilename))
        logging.info('')

        if childCount < childTrajectories:

            self.activate()
            
            # Seed the random number generator
            random_init()
    
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
        
            # Continue evolving parent trajectory, interrupting to sample sets of
            # child trajectories in order to update the recrossing factor
            parentIter = 0
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
                
                self.saveRecrossingFactor(recrossingFilename, kappa_num, kappa_denom, childCount,
                    childTrajectories, equilibrationSteps, childSamplingSteps, childEvolutionSteps, childrenPerSampling)
                
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
                kappa_num[childStep] / childCount,
            ))
        logging.info('=========== =========== ===========')
        logging.info('')

        logging.info('Final value of transmission coefficient = {0:.6f}'.format(kappa_num[-1] / kappa_denom))
        logging.info('')
        
        self.recrossingFactor = kappa_num[-1] / kappa_denom
        
        return self.recrossingFactor

    def createWorkingDirectory(self, path=None):
        """
        Create the directory used for saving the calculation output. If not
        explicitly defined, the directory is created as a subdirectory of the
        specified output directory based on the current temperature.
        """
        if path is not None:
            workingDirectory = path
        else:
            workingDirectory = os.path.join(self.outputDirectory, '{0:g}'.format(self.T))
            
        # Create the working directory on disk
        try:
            os.makedirs(workingDirectory)
        except OSError:
            pass
        
        # Return the full path to the chosen working directory 
        return os.path.abspath(workingDirectory)

    def saveUmbrellaConfigurations(self, path, evolutionSteps):
        """
        Save the results of an umbrella configurations calculation to `path` on
        disk. This serves as both a record of the calculation and a means of
        restarting an incomplete calculation.
        """
        f = open(path, 'w')
        
        f.write('****************************\n')
        f.write('RPMD umbrella configurations\n')
        f.write('****************************\n\n')

        f.write('Temperature                             = {0:g} K\n'.format(self.T))
        f.write('Number of beads                         = {0:d}\n'.format(self.Nbeads))
        f.write('Time step                               = {0:g} ps\n'.format(self.dt * 2.418884326505e-5))
        f.write('Number of umbrella integration windows  = {0:d}\n'.format(len(self.umbrellaConfigurations)))
        f.write('Trajectory evolution time               = {0:g} ps ({1:d} steps)\n\n'.format(evolutionSteps * self.dt * 2.418884326505e-5, evolutionSteps))
        
        for xi, q in self.umbrellaConfigurations:
            f.write('xi = {0:g}\n'.format(xi))
            for j in range(self.Natoms):
                f.write('{0:5} {1:11.6f} {2:11.6f} {3:11.6f}\n'.format(self.reactants.atoms[j], q[0,j], q[1,j], q[2,j]))                
            f.write('\n')
        
        f.close()
        
    def loadUmbrellaConfigurations(self, path):
        """
        Load the results of an umbrella configurations calculation from `path`
        on disk. This can be useful both as a means of postprocessing results
        at a later date and for restarting an incomplete calculation.
        """
        
        f = open(path, 'r')

        # Header
        f.readline()
        jobtype = f.readline()
        if jobtype.strip() != 'RPMD umbrella configurations':
            raise RPMDError('{0} is not a valid RPMD umbrella configurations output file.'.format(jobtype))
        f.readline()
        f.readline()
        
        # Parameters
        line = f.readline()
        while line.strip() != '':
            param, data = line.split('=')
            param = param.strip()
            data = data.split()
            if param == 'Temperature':
                T = float(data[0])
            elif param == 'Number of beads':
                Nbeads = int(data[0])
            elif param == 'Time step':
                dt = float(data[0]) / 2.418884326505e-5
            elif param == 'Number of umbrella integration windows':
                Nxi = int(data[0])
            elif param == 'Trajectory evolution time':
                evolutionSteps = int(data[2][1:])
            else:
                raise RPMDError('Invalid umbrella configurations parameter {0!r}.'.format(param))
            line = f.readline()
        
        # Configurations
        xi_list = numpy.zeros(Nxi)
        q_list = []
        for n in range(Nxi):
            line = f.readline().strip()
            xi_list[n] = float(line.split()[-1])
            q = []
            line = f.readline().strip()
            while line != '':
                label, x, y, z = line.split()
                q.append([float(x), float(y), float(z)])
                line = f.readline().strip()
            q_list.append(q)
        xi_list = numpy.array(xi_list)
        q_list = numpy.array(q_list).T
                
        f.close()
        
        return xi_list, q_list, evolutionSteps
        
    def loadUmbrellaSampling(self, path):
        """
        Load the results of an umbrella sampling calculation from `path` on
        disk. This can be useful both as a means of postprocessing results
        at a later date and for restarting an incomplete calculation.
        """
        
        f = open(path, 'r')

        # Header
        f.readline()
        jobtype = f.readline()
        if jobtype.strip() != 'RPMD umbrella sampling':
            raise RPMDError('{0} is not a valid RPMD umbrella sampling output file.'.format(jobtype))
        f.readline()
        f.readline()
        
        # Parameters
        line = f.readline()
        while line.strip() != '':
            param, data = line.split('=')
            param = param.strip()
            data = data.split()
            if param == 'Temperature':
                T = float(data[0])
            elif param == 'Number of beads':
                Nbeads = int(data[0])
            elif param == 'Time step':
                dt = float(data[0]) / 2.418884326505e-5
            elif param == 'Number of umbrella integration windows':
                Nxi = int(data[0])
            else:
                raise RPMDError('Invalid umbrella configurations parameter {0!r}.'.format(param))
            line = f.readline()
        
        # Data
        xi_list = []; av_list = []; av2_list = []; count_list = []
        line = f.readline()
        line = f.readline()
        line = f.readline()
        line = f.readline()
        while line != '' and len(line) > 8 and line[0:8] != '========':
            xi, av, av2, count, xi_mean, xi_var = line.split()
            xi_list.append(float(xi))
            av_list.append(float(av))
            av2_list.append(float(av2))
            count_list.append(int(count))
            line = f.readline().strip()
        
        xi_list = numpy.array(xi_list)
        av_list = numpy.array(av_list)
        av2_list = numpy.array(av2_list)
        count_list = numpy.array(count_list)
        
        f.close()
        
        return xi_list, av_list, av2_list, count_list
        
    def saveRecrossingFactor(self, path, kappa_num, kappa_denom, trajectoryCount,
                             childTrajectories, equilibrationSteps, 
                             childSamplingSteps, childEvolutionSteps, 
                             childrenPerSampling):
        """
        Save the results of a recrossing factor calculation to `path` on disk.
        This serves as both a record of the calculation and a means of
        restarting an incomplete calculation.
        """
        f = open(path, 'w')
        
        f.write('**********************\n')
        f.write('RPMD recrossing factor\n')
        f.write('**********************\n\n')
        
        f.write('Temperature                             = {0:g} K\n'.format(self.T))
        f.write('Number of beads                         = {0:d}\n'.format(self.Nbeads))
        f.write('Reaction coordinate                     = {0:.6f}\n'.format(self.xi_current))
        f.write('Time step                               = {0:g} ps\n'.format(self.dt * 2.418884326505e-5))
        f.write('Total number of child trajectories      = {0:d}\n'.format(childTrajectories))
        f.write('Initial parent equilibration time       = {0:g} ps ({1:d} steps)\n'.format(equilibrationSteps * self.dt * 2.418884326505e-5, equilibrationSteps))
        f.write('Frequency of child trajectory sampling  = {0:g} ps ({1:d} steps)\n'.format(childSamplingSteps * self.dt * 2.418884326505e-5, childSamplingSteps))
        f.write('Length of child trajectories            = {0:g} ps ({1:d} steps)\n'.format(childEvolutionSteps * self.dt * 2.418884326505e-5, childEvolutionSteps))
        f.write('Number of children per sampling         = {0:d}\n\n'.format(childrenPerSampling))
        
        kappa_denom = float(kappa_denom)
        
        f.write('========= ============= ============= ========= =========== ===========\n')
        f.write('Time (fs) kappa_num     kappa_denom   count     kappa (old) kappa (new)\n')
        f.write('========= ============= ============= ========= =========== ===========\n')
        for n in range(kappa_num.shape[0]):
            f.write('{0:9.3f} {1:13.4f} {2:13.4f} {3:9d} {4:11.6f} {5:11.6f}\n'.format(
                n * self.dt * 2.418884326505e-2,
                kappa_num[n],
                kappa_denom,
                trajectoryCount,
                kappa_num[n] / trajectoryCount,
                kappa_num[n] / kappa_denom,
            ))
        f.write('========= ============= ============= ========= =========== ===========\n')
        
        f.close()
    
    def loadRecrossingFactor(self, path):
        """
        Load the results of a recrossing factor calculation from `path` on disk.
        This can be useful both as a means of postprocessing results at a later
        date and for restarting an incomplete calculation.
        """

        f = open(path, 'r')

        # Header
        f.readline()
        jobtype = f.readline()
        if jobtype.strip() != 'RPMD recrossing factor':
            raise RPMDError('{0} is not a valid RPMD recrossing factor output file.')
        f.readline()
        f.readline()
        
        # Parameters
        line = f.readline()
        while line.strip() != '':
            param, data = line.split('=')
            param = param.strip()
            data = data.split()
            if param == 'Temperature':
                T = float(data[0])
            elif param == 'Number of beads':
                Nbeads = int(data[0])
            elif param == 'Reaction coordinate':
                xi_current = float(data[0])
            elif param == 'Time step':
                dt = float(data[0]) / 2.418884326505e-5
            elif param == 'Total number of child trajectories':
                childTrajectories = int(data[0])
            elif param == 'Initial parent equilibration time':
                equilibrationSteps = int(data[2][1:])
            elif param == 'Frequency of child trajectory sampling':
                childSamplingSteps = int(data[2][1:])
            elif param == 'Length of child trajectories':
                childEvolutionSteps = int(data[2][1:])
            elif param == 'Number of children per sampling':
                childrenPerSampling = int(data[0])
            else:
                raise RPMDError('Invalid recrossing factor parameter {0!r}.'.format(param))
            line = f.readline()
        
        # Data
        line = f.readline()
        line = f.readline()
        line = f.readline()
        line = f.readline().strip()
        n = 0; kappa_num = []
        while line != '' and len(line) > 8 and line[0:8] != '========':
            t, num, denom, count, kappa_old, kappa_new = line.split()
            kappa_num.append(float(num))
            kappa_denom = float(denom)
            trajectoryCount = int(count)
            line = f.readline().strip()
        kappa_num = numpy.array(kappa_num, order='F')
        kappa_denom = numpy.array(kappa_denom, order='F')
        
        f.close()

        return (T, Nbeads, xi_current, dt, kappa_num, kappa_denom, trajectoryCount,
            childTrajectories, equilibrationSteps, childSamplingSteps, 
            childEvolutionSteps, childrenPerSampling)
    
    def computeRPMDRateCoefficient(self):
        """
        Compute the value of the RPMD rate coefficient.
        """
        
        if self.potentialOfMeanForce is None or self.recrossingFactor is None:
            raise RPMDError('You must run computePotentialOfMeanForce() and computeRecrossingFactor() before running computeRPMDRateCoefficient().')
        
        logging.info('*********************')
        logging.info('RPMD rate coefficient')
        logging.info('*********************')
        logging.info('')

        logging.info('Parameters')
        logging.info('==========')
        logging.info('Temperature                             = {0:g} K'.format(self.T))
        logging.info('')
        
        fromAtomicUnits = 1e6 * ((5.2917721092e-11)**3 / 2.418884326505e-17) * constants.Na
        
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
