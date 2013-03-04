#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RPMDrate - Bimolecular reaction rates via ring polymer molecular dynamics
#
#   Copyright (c) 2012 by Joshua W. Allen (jwallen@mit.edu)
#                         William H. Green (whgreen@mit.edu)
#                         Yury V. Suleimanov (ysuleyma@mit.edu, ysuleyma@princeton.edu)
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

import rpmdrate.constants as constants
import rpmdrate.quantity as quantity

from rpmdrate._main import *
from rpmdrate.surface import TransitionState

################################################################################

class RPMDError(Exception):
    """
    An exception raised when an error occurs during an RPMD simulation. Pass a
    string describing the circumstances of the exceptional behavior.
    """
    pass

################################################################################

def runUmbrellaTrajectory(rpmd, xi_current, p, q, equilibrationSteps, evolutionSteps, kforce, xi_range, saveTrajectory):
    """
    Run an individual umbrella integration trajectory, returning the sum of the
    first and second moments of the reaction coordinate at each time step.
    """
    rpmd.activate()
    steps = 0
    while steps < evolutionSteps:
        p1 = numpy.asfortranarray(p.copy())
        q1 = numpy.asfortranarray(q.copy())
        result = system.equilibrate(0, p1, q1, equilibrationSteps, xi_current, rpmd.potential, kforce, False, saveTrajectory)
        if result != 0: continue
        dav, dav2, actualSteps, result = system.umbrella_trajectory(0, p1, q1, evolutionSteps - steps, xi_current, rpmd.potential, kforce, xi_range, saveTrajectory)
        steps += actualSteps
        if result != 0: continue
    
    return dav, dav2, steps, p1, q1

def runRecrossingTrajectory(rpmd, xi_current, p, q, evolutionSteps, saveTrajectory):
    """
    Run an individual pair of recrossing factor child trajectories, returning
    the contributions to the numerator and denominator of the recrossing factor
    from this trajectory pair. We use pairs of trajectories so that we always
    sample in the positive and negative directions of the initial sampled
    momenta.
    """
    rpmd.activate()
    result1 = 1; result2 = 1
    while result1 != 0 or result2 != 0:
        # Trajectory for the negative of the sampled momenta
        t1 = numpy.array(0.0, order='F')
        p1 = (-p).copy('F')
        q1 = q.copy('F')
        kappa_num1 = numpy.zeros(evolutionSteps, order='F')
        kappa_denom1 = numpy.array(0.0, order='F')
        result1 = system.recrossing_trajectory(t1, p1, q1, xi_current, rpmd.potential, saveTrajectory, kappa_num1, kappa_denom1)
        if result1 != 0: continue

        # Trajectory for the positive of the sampled momenta
        t2 = numpy.array(0.0, order='F')
        p2 = p.copy('F')
        q2 = q.copy('F')
        kappa_num2 = numpy.zeros(evolutionSteps, order='F')
        kappa_denom2 = numpy.array(0.0, order='F')
        result2 = system.recrossing_trajectory(t2, p2, q2, xi_current, rpmd.potential, saveTrajectory, kappa_num2, kappa_denom2)
        if result2 != 0: continue
    
    return kappa_num1 + kappa_num2, kappa_denom1 + kappa_denom2

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
    `xi_range`                  The allowed maximum valid range of the reaction coordinate in this window
    --------------------------- ------------------------------------------------
    `count`                     The number of samples taken
    `av`                        The mean of the reaction coordinate times the number of samples
    `av2`                       The variance of the reaction coordinate times the number of samples
    =========================== ================================================    
    
    """
    
    def __init__(self, xi=None, kforce=None, trajectories=None, equilibrationTime=None, evolutionTime=None, xi_range=0.0):
        # These parameters control the umbrella sampling trajectories
        self.xi = xi
        self.kforce = kforce
        self.trajectories = trajectories
        if equilibrationTime is not None:
            self.equilibrationTime = float(quantity.convertTime(equilibrationTime, "ps")) / 2.418884326505e-5
        else:
            self.equilibrationTime = None
        if evolutionTime is not None:
            self.evolutionTime = float(quantity.convertTime(evolutionTime, "ps")) / 2.418884326505e-5
        else:
            self.evolutionTime = None
        self.xi_range = xi_range
        self.q = None
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

    def __init__(self, label, T, Nbeads, reactants, transitionState, potential, thermostat, processes=1, outputDirectory='.', randomSeed=None):
        """
        Initialize an RPMD object. The `mass` of each atom should be given in
        g/mol, while the `Rinf` value should be given in angstroms. (They will
        be converted to atomic units.)
        """
        self.label = label
        self.T = T
        self.Nbeads = Nbeads
        self.mass = reactants.mass
        self.Natoms = len(self.mass)
        self.reactants = reactants
        self.transitionStates = [transitionState]
        self.potential = potential
        self.thermostat = thermostat
        self.processes = processes
        self.outputDirectory = os.path.abspath(outputDirectory)
        self.randomSeed = randomSeed
        
        self.beta = 4.35974417e-18 / (constants.kB * self.T)
        self.dt = 0
        self.xi_current = 0
        self.mode = 0
        
        self.umbrellaConfigurations = None
        self.umbrellaWindows = None
        self.potentialOfMeanForce = None
        self.recrossingFactor = None
    
    def addEquivalentTransitionState(self, formingBonds, breakingBonds):
        """
        Add an equivalent transition state to the RPMD system, defined by lists
        of forming and breaking bonds `formingBonds` and `breakingBonds`,
        respectively. The bonds must correspond to the forming and breaking
        bonds in the original transition state.
        """
        mapping = {}
        for bond1, bond2 in zip(self.transitionStates[0].formingBonds, formingBonds):
            atom11, atom12 = bond1
            atom21, atom22 = bond2
            if atom11 in mapping and mapping[atom11] != atom21:
                raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom11, mapping[atom11], atom21))
            elif atom21 in mapping and mapping[atom21] != atom11:
                raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom21, mapping[atom21], atom11))
            else:
                mapping[atom11] = atom21
                mapping[atom21] = atom11
            if atom12 in mapping and mapping[atom12] != atom22:
                raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom12, mapping[atom12], atom22))
            elif atom22 in mapping and mapping[atom22] != atom12:
                raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom22, mapping[atom22], atom12))
            else:
                mapping[atom12] = atom22
                mapping[atom22] = atom12
        for bond1, bond2 in zip(self.transitionStates[0].breakingBonds, breakingBonds):
            atom11, atom12 = bond1
            atom21, atom22 = bond2
            if atom11 in mapping and mapping[atom11] != atom21:
                raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom11, mapping[atom11], atom21))
            elif atom21 in mapping and mapping[atom21] != atom11:
                raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom21, mapping[atom21], atom11))
            else:
                mapping[atom11] = atom21
                mapping[atom21] = atom11
            if atom12 in mapping and mapping[atom12] != atom22:
                raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom12, mapping[atom12], atom22))
            elif atom22 in mapping and mapping[atom22] != atom12:
                raise ValueError('Inconsistent indices in equivalent transition state: {0} mapped to both {1} and {2}.'.format(atom22, mapping[atom22], atom12))
            else:
                mapping[atom12] = atom22
                mapping[atom22] = atom12
        for atom in range(1, self.Natoms+1):
            if atom not in mapping:
                mapping[atom] = atom
        
        geometry = numpy.zeros_like(self.transitionStates[0].geometry)
        for atom in range(self.Natoms):
            geometry[:,atom] = self.transitionStates[0].geometry[:,mapping[atom+1]-1]
        
        self.transitionStates.append(TransitionState(
            geometry = (geometry.T,"bohr"),
            formingBonds = formingBonds,
            breakingBonds = breakingBonds,
        ))
    
    def activate(self, Nbeads=None):
        """
        Set this object as the active RPMD system in the Fortran layer. Note
        that the dividing surface properties must be set in the corresponding
        modules in ``rpmd._main``, *not* in ``rpmd._surface``.
        """
        Natoms = self.mass.shape[0]
        Nbeads = self.Nbeads if Nbeads is None else Nbeads
        system.dt = self.dt
        system.beta = self.beta
        system.mass[0:Natoms] = self.mass
        system.mode = self.mode
        self.reactants.activate(module=reactants)
        self.thermostat.activate(module=system, Natoms=Natoms, Nbeads=Nbeads)
        
        Nts = len(self.transitionStates)
        Nforming_bonds = max([ts.formingBonds.shape[0] for ts in self.transitionStates])
        Nbreaking_bonds = max([ts.breakingBonds.shape[0] for ts in self.transitionStates])
        assert Nforming_bonds == Nbreaking_bonds
        
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
        transition_state.number_of_bonds = Nforming_bonds
        transition_state.forming_bonds[0:Nts,0:Nforming_bonds,:] = formingBonds
        transition_state.forming_bond_lengths[0:Nts,0:Nforming_bonds] = formingBondLengths
        transition_state.breaking_bonds[0:Nts,0:Nbreaking_bonds,:] = breakingBonds
        transition_state.breaking_bond_lengths[0:Nts,0:Nbreaking_bonds] = breakingBondLengths
    
    def generateUmbrellaConfigurations(self, 
                                       dt, 
                                       evolutionTime,
                                       xi_list,
                                       kforce):
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
        Nbeads = 1
        xi_list = numpy.array(xi_list)
        Nxi = len(xi_list)
        geometry = self.transitionStates[0].geometry
        thermostat = self.thermostat
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
        logging.info('Number of beads                         = {0:d}'.format(Nbeads))
        logging.info('Time step                               = {0:g} ps'.format(self.dt * 2.418884326505e-5))
        logging.info('Number of umbrella integration windows  = {0:d}'.format(Nxi))
        logging.info('Trajectory evolution time               = {0:g} ps ({1:d} steps)'.format(evolutionSteps * self.dt * 2.418884326505e-5, evolutionSteps))
        logging.info('')

        # Set up output files and directory
        workingDirectory = self.createWorkingDirectory()
        # The umbrella configurations are independent of temperature and number
        # of beads, so store them in the top-level directory
        configurationsFilename = os.path.realpath(os.path.join(workingDirectory, '..', '..', 'umbrella_configurations.dat'))

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
        self.activate(Nbeads)

        # Seed the random number generator
        self.initializeRandomNumberGenerator()

        # Generate initial position using transition state geometry
        # (All beads start at same position)
        q = numpy.zeros((3,self.Natoms,Nbeads), order='F')
        for i in range(3):
            for j in range(self.Natoms):
                for k in range(Nbeads):
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
            logging.info('Generating configuration at xi = {0:.4f} for {1:g} ps...'.format(xi_current, evolutionSteps * self.dt * 2.418884326505e-5))
            p = self.sampleMomentum(Nbeads=Nbeads)
            result = system.equilibrate(0, p, q, evolutionSteps, xi_current, self.potential, kforce[l], False, False)
            logging.info('Finished generating configuration at xi = {0:.4f}.'.format(xi_current))
            q_initial[:,:,l] = q[:,:,0]
                        
        # Now start at xi = 1 and move in the xi < 1 direction, using the
        # result of the previous xi as the initial position for the next xi
        q[:,:,0] = q_initial[:,:,start]
        for l in range(start - 1, -1, -1):
            xi_current = xi_list[l]
            
            # Equilibrate in this window
            logging.info('Generating configuration at xi = {0:.4f} for {1:g} ps...'.format(xi_current, evolutionSteps * self.dt * 2.418884326505e-5))
            p = self.sampleMomentum(Nbeads=Nbeads)
            result = system.equilibrate(0, p, q, evolutionSteps, xi_current, self.potential, kforce[l], False, False)
            logging.info('Finished generating configuration at xi = {0:.4f}.'.format(xi_current))
            q_initial[:,:,l] = q[:,:,0]
        
        # Store the computed configurations on the object for future use in
        # umbrella sampling
        self.umbrellaConfigurations = []
        for l in range(Nxi):
            xi_current = xi_list[l]
            q_current = self.cleanGeometry(q_initial[:,:,l])
            self.umbrellaConfigurations.append((xi_current, q_current))
            logging.info('Configuration at xi = {0:.4f}:'.format(xi_current))
            for j in range(self.Natoms):
                logging.info('{0:5} {1:11.6f} {2:11.6f} {3:11.6f}'.format(self.reactants.atoms[j], q_current[0,j], q_current[1,j], q_current[2,j]))                
            logging.info('')
        
        self.saveUmbrellaConfigurations(configurationsFilename, evolutionSteps)
    
    def conductUmbrellaSampling(self, 
                                dt, 
                                windows,
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
        Nwindows = len(windows)
        thermostat = self.thermostat
        self.mode = 1
        processes = self.processes
        
        # Create a pool of subprocesses to farm out the individual trajectories to
        try:
            import multiprocessing
            if processes > 1:
                pool = multiprocessing.Pool(processes=processes)
            else:
                pool = None
        except ImportError:
            if processes > 1:
                raise ValueError('The "multiprocessing" package was not found in this Python installation; you must install this package or set processes to 1.')
            pool = None
        results = []

        logging.info('**********************')
        logging.info('RPMD umbrella sampling')
        logging.info('**********************')
        logging.info('')
        
        logging.info('Parameters')
        logging.info('==========')
        logging.info('Temperature                             = {0:g} K'.format(self.T))
        logging.info('Number of beads                         = {0:d}'.format(self.Nbeads))
        logging.info('Time step                               = {0:g} ps'.format(self.dt * 2.418884326505e-5))
        logging.info('Number of umbrella integration windows  = {0:d}'.format(Nwindows))
        logging.info('')

        # Set up output files and directory
        workingDirectory = self.createWorkingDirectory()

        self.activate()

        # Seed the random number generator
        self.initializeRandomNumberGenerator()

        # Load any previous umbrella sampling trajectories for each window
        for window in windows:
            umbrellaFilename = os.path.join(workingDirectory, 'umbrella_sampling_{0:.4f}.dat'.format(window.xi))
            if os.path.exists(umbrellaFilename):
                # Previous trajectories existed, so load them
                logging.info('Loading saved output for xi = {0:.4f} from {1}'.format(window.xi, umbrellaFilename))
                xi, kforce, av_list, av2_list, count_list = self.loadUmbrellaSampling(umbrellaFilename)
                assert abs(xi - window.xi) < 1e-6
                if len(av_list) > 0:
                    window.av += av_list[-1]
                    window.av2 += av2_list[-1]
                    window.count += count_list[-1]
                    
            else:
                # No previous trajectories existed, so start a new output file
                # for this window
                equilibrationSteps = int(round(window.equilibrationTime / self.dt))
                evolutionSteps = int(round(window.evolutionTime / self.dt))
                f = open(umbrellaFilename, 'w')
                f.write('**********************\n')
                f.write('RPMD umbrella sampling\n')
                f.write('**********************\n\n')
                f.write('Temperature                             = {0:g} K\n'.format(self.T))
                f.write('Number of beads                         = {0:d}\n'.format(self.Nbeads))
                f.write('Time step                               = {0:g} ps\n'.format(self.dt * 2.418884326505e-5))
                f.write('Reaction coordinate                     = {0:.4f}\n'.format(window.xi))
                f.write('Equilibration time                      = {0:g} ps ({1:d} steps)\n'.format(equilibrationSteps * self.dt * 2.418884326505e-5, equilibrationSteps))
                f.write('Trajectory evolution time               = {0:g} ps ({1:d} steps)\n'.format(evolutionSteps * self.dt * 2.418884326505e-5, evolutionSteps))
                f.write('Force constant                          = {0:g}\n'.format(window.kforce))
                if window.xi_range:
                    f.write('Valid reaction coordinate range         = {0:g}\n'.format(window.xi_range))
                f.write('\n')
                    
                f.write('=============== =============== =========== =============== ===============\n')
                f.write('total av        total av2       count       xi_mean         xi_var\n')
                f.write('=============== =============== =========== =============== ===============\n')
                f.flush()
                os.fsync(f.fileno())
                f.close()
                
        # Load initial configuration using results from generateUmbrellaConfigurations()
        for window in windows:
            window.q = numpy.empty((3,self.Natoms,self.Nbeads), order='F')
            for xi, q_initial in self.umbrellaConfigurations:
                if xi >= window.xi:
                    break
            for k in range(self.Nbeads):
                window.q[:,:,k] = q_initial

        # This implementation is breadth-first, as we would rather get some
        # data in all windows than get lots of data in a few windows
        done = False
        while not done:
            
            # Clear trajectories from previous iteration
            results = []
            done = True
            
            # Run one trajectory for each window that needs more sampling
            for window in windows:

                equilibrationSteps = int(round(window.equilibrationTime / self.dt))
                evolutionSteps = int(round(window.evolutionTime / self.dt))
            
                if window.count >= window.trajectories * evolutionSteps:
                    # We've already sampled enough in this window, so don't
                    # do any more
                    continue
                
                done = False
                
                # Load initial configuration from last finished trajectory
                # for this window
                q = window.q[:,:,:]
                
                # Spawn sampling trajectory in this window
                windowEvolutionSteps = evolutionSteps
                windowEquilibrationSteps = equilibrationSteps
                logging.info('Spawning sampling trajectory at xi = {0:.4f}...'.format(window.xi))
                p = self.sampleMomentum()
                args = (self, window.xi, p, q, windowEquilibrationSteps, windowEvolutionSteps, window.kforce, window.xi_range, saveTrajectories)
                if pool:
                    results.append([window, pool.apply_async(runUmbrellaTrajectory, args)])
                else:
                    results.append([window, runUmbrellaTrajectory(*args)])
              
            count = 0  
            for window, result in results:

                logging.info('Processing trajectory at xi = {0:.4f}...'.format(window.xi))
                    
                # This line will block until the trajectory finishes
                if pool:
                    dav, dav2, dcount, p, q = result.get()
                else:
                    dav, dav2, dcount, p, q = result
                
                # Update the mean and variance with the results from this trajectory
                # Note that these are counted at each time step in each trajectory
                window.av += dav
                window.av2 += dav2
                window.count += dcount
                
                # Print the updated mean and variance to the log file
                av = window.av / window.count
                av2 = window.av2 / window.count
                mean = av
                variance = av2 - av * av
                if dcount > 0:
                    logging.info('{0:11d} {1:15.8f} {2:15.8f} {3:15.5e}'.format(window.count, av, av2, variance))
    
                # Also save the geometry to use as the initial geometry for
                # the next trajectory in this window
                window.q = q[:,:,:]
    
                umbrellaFilename = os.path.join(workingDirectory, 'umbrella_sampling_{0:.4f}.dat'.format(window.xi))
                f = open(umbrellaFilename, 'a')
                f.write('{0:15.8f} {1:15.8f} {2:11d} {3:15.5e} {4:15.5e}\n'.format(window.av, window.av2, window.count, mean, variance))
                f.flush()
                os.fsync(f.fileno())
                f.close()
                
                count += 1
                        
        logging.info('')
        
    def computePotentialOfMeanForce(self, windows=None, xi_min=None, xi_max=None, bins=5000):
        """
        Compute the potential of mean force of the system at the given
        temperature by integrating over the given reaction coordinate range
        using the given number of bins. This requires that you have previously
        used umbrella sampling to determine the mean and variance in each bin.
        """               
        # Set up output files and directory
        workingDirectory = self.createWorkingDirectory()
        potentialFilename = os.path.join(workingDirectory, 'potential_of_mean_force.dat')

        # If windows is not specified, then try to determine the available
        # windows by loading from the files in the working directory
        if windows is None:
            windows = []
            for root, dirs, files in os.walk(workingDirectory):
                for f in files:
                    if f.startswith('umbrella_sampling_'):
                        umbrellaFilename = os.path.join(root, f)
                        logging.info('Loading saved output from {0}'.format(umbrellaFilename))
                        xi, kforce, av_list, av2_list, count_list = self.loadUmbrellaSampling(umbrellaFilename)
                        if len(av_list) > 0:
                            window = Window(xi=xi, kforce=kforce)
                            window.av += av_list[-1]
                            window.av2 += av2_list[-1]
                            window.count += count_list[-1]
                            windows.append(window)
            if len(windows) == 0:
                raise RPMDError('No windows specified to computePotentialOfMeanForce(), and could not determine the available windows from the saved files in the working directory.')
            windows.sort(key=lambda w: w.xi)
            
        # If xi_min and/or xi_max are not given, then default to the entire
        # range of windows
        if xi_min is None: xi_min = windows[0].xi
        if xi_max is None: xi_max = windows[-1].xi
            
        self.umbrellaWindows = windows
        Nwindows = len(self.umbrellaWindows)
        
        xi_list = numpy.linspace(xi_min, xi_max, bins, True)
        
        logging.info('****************************')
        logging.info('RPMD potential of mean force')
        logging.info('****************************')
        logging.info('')
        
        logging.info('Parameters')
        logging.info('==========')
        logging.info('Number of umbrella integration windows  = {0:d}'.format(Nwindows))
        logging.info('Lower bound of reaction coordinate      = {0:g}'.format(xi_min))
        logging.info('Upper bound of reaction coordinate      = {0:g}'.format(xi_max))
        logging.info('Number of bins                          = {0:d}'.format(bins))
        logging.info('')

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
        self.potentialOfMeanForce = numpy.zeros((2,bins-1))
        A = 0.0
        for n in range(bins-1):
            dx = xi_list[n+1] - xi_list[n]
            self.potentialOfMeanForce[0,n] = 0.5 * (xi_list[n] + xi_list[n+1])
            A += 0.5 * dx * (dA[n] + dA[n+1])
            self.potentialOfMeanForce[1,n] = A
        self.potentialOfMeanForce[1,:] -= numpy.min(self.potentialOfMeanForce[1,:])
        
        # Save the results to file
        self.savePotentialOfMeanForce(potentialFilename)

    def computeRecrossingFactor(self, 
                                dt, 
                                equilibrationTime,
                                childTrajectories,
                                childSamplingTime,
                                childrenPerSampling,
                                childEvolutionTime,
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
        self.kforce = 0.0
        self.xi_current = xi_current
        thermostat = self.thermostat
        self.mode = 2
        processes = self.processes
        
        # Load the geometry from the closest umbrella configuration
        for xi, q in self.umbrellaConfigurations:
            if xi >= xi_current:
                geometry = q
                break
        else:
            geometry = self.transitionStates[0].geometry
        
        # Initialize parameters used to compute recrossing factor
        kappa_num = numpy.zeros(childEvolutionSteps, order='F')
        kappa_denom = numpy.array(0.0, order='F')
        childCount = 0
        
        # Create a pool of subprocesses to farm out the individual trajectories to
        try:
            import multiprocessing
            if processes > 1:
                pool = multiprocessing.Pool(processes=processes)
            else:
                pool = None
        except ImportError:
            if processes > 1:
                raise ValueError('The "multiprocessing" package was not found in this Python installation; you must install this package or set processes to 1.')
            pool = None
        results = []

        logging.info('**********************')
        logging.info('RPMD recrossing factor')
        logging.info('**********************')
        logging.info('')
        
        logging.info('Parameters')
        logging.info('==========')
        logging.info('Temperature                             = {0:g} K'.format(self.T))
        logging.info('Number of beads                         = {0:d}'.format(self.Nbeads))
        logging.info('Reaction coordinate                     = {0:.4f}'.format(xi_current))
        logging.info('Time step                               = {0:g} ps'.format(self.dt * 2.418884326505e-5))
        logging.info('Total number of child trajectories      = {0:d}'.format(childTrajectories))
        logging.info('Initial parent equilibration time       = {0:g} ps ({1:d} steps)'.format(equilibrationSteps * self.dt * 2.418884326505e-5, equilibrationSteps))
        logging.info('Frequency of child trajectory sampling  = {0:g} ps ({1:d} steps)'.format(childSamplingSteps * self.dt * 2.418884326505e-5, childSamplingSteps))
        logging.info('Length of child trajectories            = {0:g} ps ({1:d} steps)'.format(childEvolutionSteps * self.dt * 2.418884326505e-5, childEvolutionSteps))
        logging.info('Number of children per sampling         = {0:d}'.format(childrenPerSampling))
        logging.info('')
        
        # Set up output files and directory
        workingDirectory = self.createWorkingDirectory()
        recrossingFilename = os.path.join(workingDirectory, 'recrossing_factor_{0:.4f}.dat'.format(self.xi_current))

        # Look for existing output file for this calculation
        # If a file exists, we won't repeat the calculation unless more
        # child trajectories are requested
        if os.path.exists(recrossingFilename):
            logging.info('Loading saved output from {0}'.format(recrossingFilename))
            (T0, Nbeads0, xi_current0, dt0, kappa_num0, kappa_denom0, trajectoryCount0,
                childTrajectories0, equilibrationSteps0, childSamplingSteps0, 
                childEvolutionSteps0, childrenPerSampling0) = self.loadRecrossingFactor(recrossingFilename)
            if T0 == self.T and Nbeads0 == self.Nbeads and dt0 == self.dt and abs(xi_current0 - self.xi_current) < 1e-4:
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

        recrossingFactor = []
        
        if childCount < childTrajectories:

            self.activate()
            
            # Seed the random number generator
            self.initializeRandomNumberGenerator()

            # Generate initial position using transition state geometry
            # (All beads start at same position)
            q0 = numpy.zeros((3,self.Natoms,self.Nbeads), order='F')
            for i in range(3):
                for j in range(self.Natoms):
                    for k in range(self.Nbeads):
                        q0[i,j,k] = geometry[i,j]
            
            # Equilibrate parent trajectory while constraining to dividing surface
            # and sampling from Andersen thermostat
            logging.info('Equilibrating parent trajectory for {0:g} ps...'.format(equilibrationSteps * self.dt * 2.418884326505e-5))
            result = 1
            while result != 0:
                q = numpy.asfortranarray(q0.copy())
                p = self.sampleMomentum()            
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
                    if pool:
                        results.append(pool.apply_async(runRecrossingTrajectory, args))
                    else:
                        results.append(runRecrossingTrajectory(*args))           
                    childCount += 2
    
                for child in range(childrenPerSampling / 2):
                    # This line will block until the child trajectory finishes
                    if pool:
                        num, denom = results[child].get()
                    else:
                        num, denom = results[child]
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
            
            logging.info('Finished sampling of {0:d} child trajectories.'.format(childCount))
            logging.info('')
        
        logging.info('Result of recrossing factor calculation:')
        logging.info('')
        logging.info('=========== ===========')
        logging.info('Time (fs)   kappa')
        logging.info('=========== ===========')
        
        for childStep in range(childEvolutionSteps):
            logging.info('{0:11.3f} {1:11.6f}'.format(
                childStep * self.dt * 2.418884326505e-2,
                kappa_num[childStep] / kappa_denom,
            ))
        logging.info('=========== ===========')
        logging.info('')

        logging.info('Final value of recrossing factor = {0:.6f}'.format(kappa_num[-1] / kappa_denom))
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
            workingDirectory = os.path.join(self.outputDirectory, '{0:g}'.format(self.T), '{0:d}'.format(self.Nbeads))
            
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
        f.write('Number of beads                         = 1\n')
        f.write('Time step                               = {0:g} ps\n'.format(self.dt * 2.418884326505e-5))
        f.write('Number of umbrella integration windows  = {0:d}\n'.format(len(self.umbrellaConfigurations)))
        f.write('Trajectory evolution time               = {0:g} ps ({1:d} steps)\n\n'.format(evolutionSteps * self.dt * 2.418884326505e-5, evolutionSteps))
        
        for xi, q in self.umbrellaConfigurations:
            f.write('xi = {0:.4f}\n'.format(xi))
            for j in range(self.Natoms):
                f.write('{0:5} {1:11.6f} {2:11.6f} {3:11.6f}\n'.format(self.reactants.atoms[j], q[0,j], q[1,j], q[2,j]))                
            f.write('\n')
        
        f.flush()
        os.fsync(f.fileno())

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
            elif param == 'Reaction coordinate':
                xi = float(data[0])
            elif param == 'Equilibration time':
                equilibrationSteps = int(data[2][1:])
            elif param == 'Trajectory evolution time':
                evolutionSteps = int(data[2][1:])
            elif param == 'Force constant':
                kforce = float(data[0])
            elif param == 'Valid reaction coordinate range':
                xi_range = float(data[0])
            else:
                raise RPMDError('Invalid umbrella sampling parameter {0!r}.'.format(param))
            line = f.readline()
        
        # Data
        av_list = []; av2_list = []; count_list = []
        line = f.readline()
        line = f.readline()
        line = f.readline()
        line = f.readline()
        while line != '' and len(line) > 8 and line[0:8] != '========':
            av, av2, count, xi_mean, xi_var = line.split()
            av_list.append(float(av))
            av2_list.append(float(av2))
            count_list.append(int(count))
            line = f.readline().strip()
        
        av_list = numpy.array(av_list)
        av2_list = numpy.array(av2_list)
        count_list = numpy.array(count_list)
        
        f.close()
        
        return xi, kforce, av_list, av2_list, count_list
        
    def savePotentialOfMeanForce(self, path):
        """
        Save the results of a potential of mean force calculation to `path` on
        disk. This serves as both a record of the calculation and a means of
        restarting an incomplete calculation.
        """
        
        xi_min = self.potentialOfMeanForce[0,0]
        xi_max = self.potentialOfMeanForce[0,-1]
        bins = self.potentialOfMeanForce.shape[1]
        
        f = open(path, 'w')
        
        f.write('****************************\n')
        f.write('RPMD potential of mean force\n')
        f.write('****************************\n\n')
        
        f.write('Temperature                             = {0:g} K\n'.format(self.T))
        f.write('Lower bound of reaction coordinate      = {0:g}\n'.format(xi_min))
        f.write('Upper bound of reaction coordinate      = {0:g}\n'.format(xi_max))
        f.write('Number of bins                          = {0:d}\n\n'.format(bins))

        f.write('=========== ===============\n')
        f.write('Rxn coord   PMF (eV)\n')
        f.write('=========== ===============\n')
        for n in range(0, self.potentialOfMeanForce.shape[1]):
            f.write('{0:11.6f} {1:11.10f}\n'.format(
                self.potentialOfMeanForce[0,n],
                self.potentialOfMeanForce[1,n] * 27.211,
            ))
        f.write('=========== ===============\n')

        f.flush()
        os.fsync(f.fileno())

        f.close()
    
    def loadPotentialOfMeanForce(self, path):
        """
        Load the results of a potential of mean force calculation from `path` on
        disk. This can be useful both as a means of postprocessing results at a
        later date and for restarting an incomplete calculation.
        """
        
        f = open(path, 'r')

        # Header
        f.readline()
        jobtype = f.readline()
        if jobtype.strip() != 'RPMD potential of mean force':
            raise RPMDError('{0} is not a valid RPMD potential of mean force output file.')
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
            elif param == 'Lower bound of reaction coordinate':
                xi_min = float(data[0])
            elif param == 'Upper bound of reaction coordinate':
                xi_max = float(data[0])
            elif param == 'Number of bins':
                bins = int(data[0])
            else:
                raise RPMDError('Invalid potential of mean force parameter {0!r}.'.format(param))
            line = f.readline()
        
        # Data
        self.potentialOfMeanForce = numpy.zeros((2,bins))
        line = f.readline()
        line = f.readline()
        line = f.readline()
        for n in range(bins):
            xi, A = f.readline().strip().split()
            self.potentialOfMeanForce[0,n] = float(xi)
            self.potentialOfMeanForce[1,n] = float(A) / 27.211
        line = f.readline()

        f.close()
    
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
        f.write('Reaction coordinate                     = {0:.4f}\n'.format(self.xi_current))
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
        
        f.flush()
        os.fsync(f.fileno())

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
    
    def saveRateCoefficient(self, path, k_QTST_s0, staticFactor, k_QTST, recrossingFactor, k_RPMD):
        """
        Save the results of a rate coefficient calculation to `path` on disk.
        This serves as both a record of the calculation and a means of
        restarting an incomplete calculation.
        """
        
        fromAtomicUnits = 1e6 * ((5.2917721092e-11)**3 / 2.418884326505e-17)

        f = open(path, 'w')
        
        f.write('*********************\n')
        f.write('RPMD rate coefficient\n')
        f.write('*********************\n\n')
        logging.info('')

        f.write('Temperature                             = {0:g} K\n'.format(self.T))
        f.write('Number of beads                         = {0:d}\n\n'.format(self.Nbeads))

        f.write('k_QTST(T;s0)                            = {0:g} cm^3/(molecule*s)\n'.format(k_QTST_s0 * fromAtomicUnits))        
        f.write('                                        = {0:g} cm^3/(mol*s)\n\n'.format(k_QTST_s0 * fromAtomicUnits * constants.Na))        
        
        f.write('Static factor                           = {0:g}\n\n'.format(staticFactor))
        
        f.write('Maximum reaction coordinate (xi_max)    = {0:.6f}\n\n'.format(self.xi_current))

        f.write('k_QTST(T;xi_max)                        = {0:g} cm^3/(molecule*s)\n'.format(k_QTST * fromAtomicUnits))
        f.write('                                        = {0:g} cm^3/(mol*s)\n\n'.format(k_QTST * fromAtomicUnits * constants.Na))
        
        f.write('Recrossing factor                       = {0:g}\n\n'.format(recrossingFactor))
        
        f.write('k_RPMD(T)                               = {0:g} cm^3/(molecule*s)\n'.format(k_RPMD * fromAtomicUnits))
        f.write('                                        = {0:g} cm^3/(mol*s)\n\n'.format(k_RPMD * fromAtomicUnits * constants.Na))
        
        f.flush()
        os.fsync(f.fileno())

        f.close()

    def loadRateCoefficient(self, path):
        """
        Load the results of a rate coefficient calculation from `path` on disk.
        This can be useful both as a means of postprocessing results at a later
        date and for restarting an incomplete calculation.
        """

        f = open(path, 'r')

        toAtomicUnits = 1e-6 / ((5.2917721092e-11)**3 / 2.418884326505e-17)

        # Header
        f.readline()
        jobtype = f.readline()
        if jobtype.strip() != 'RPMD rate coefficient':
            raise RPMDError('{0} is not a valid RPMD rate coefficient output file.')
        f.readline()
        f.readline()
        
        # Parameters
        line = f.readline()
        while line != '':
            param, data = line.split('=')
            param = param.strip()
            data = data.split()
            if param == 'Temperature':
                T = float(data[0])
            elif param == 'Number of beads':
                Nbeads = int(data[0])
            elif param == 'k_QTST(T;s0)':
                k_QTST_s0 = float(data[0]) * toAtomicUnits
            elif param == 'Static factor':
                staticFactor = float(data[0])
            elif param == 'Maximum reaction coordinate (xi_max)':
                xi_current = float(data[0])
            elif param == 'k_QTST(T;xi_max)':
                k_QTST = float(data[0]) * toAtomicUnits
            elif param == 'Recrossing factor':
                recrossingFactor = float(data[0])
            elif param == 'k_RPMD(T)':
                k_RPMD = float(data[0]) * toAtomicUnits
            elif param != '':
                raise RPMDError('Invalid recrossing factor parameter {0!r}.'.format(param))
            line = f.readline()
        
        f.close()

        return (T, Nbeads, k_QTST_s0, staticFactor, xi_current, k_QTST, recrossingFactor, k_RPMD)

    def computeRateCoefficient(self):
        """
        Compute the value of the RPMD rate coefficient.
        """
        
        if self.potentialOfMeanForce is None:
            raise RPMDError('You must run computePotentialOfMeanForce() before running computeRateCoefficient().')
        
        logging.info('*********************')
        logging.info('RPMD rate coefficient')
        logging.info('*********************')
        logging.info('')

        logging.info('Parameters')
        logging.info('==========')
        logging.info('Temperature                             = {0:g} K'.format(self.T))
        logging.info('Number of beads                         = {0:d}'.format(self.Nbeads))
        logging.info('')
        
        if self.recrossingFactor is None:
            logging.warning('computeRecrossingFactor() not before running computeRateCoefficient(); recrossing factor assumed to be unity.')
            index = numpy.argmax(self.potentialOfMeanForce[1,:])
            xi_current = self.potentialOfMeanForce[0,index]
        else:
            xi_current = self.xi_current
        
        # Set up output files and directory
        workingDirectory = self.createWorkingDirectory()
        rateFilename = os.path.join(workingDirectory, 'rate_coefficient_{0:.4f}.dat'.format(self.xi_current))

        # Compute the rate coefficient using the reactant dividing surface
        Rinf = self.reactants.Rinf
        mA = self.reactants.totalMass1
        mB = self.reactants.totalMass2
        mu = mA * mB / (mA + mB)
        k_QTST_s0 = float(4 * constants.pi * Rinf * Rinf / numpy.sqrt(2 * constants.pi * self.beta * mu))
        
        # Compute the static factor
        # (Use the same xi_current as used in the recrossing factor calculation)
        from rpmdrate.interpolate import LinearInterpolator
        f = LinearInterpolator(self.potentialOfMeanForce[0,:], self.potentialOfMeanForce[1,:])
        W1 = f(xi_current)
        W0 = f(0.0)
        staticFactor = float(numpy.exp(-self.beta * (W1 - W0)))
        
        # Correct the rate coefficient to the transition state dividing surface
        k_QTST = k_QTST_s0 * staticFactor
        
        # Correct the rate coefficient for recrossings
        if self.recrossingFactor is None:
            recrossingFactor = 1.0
        else:
            recrossingFactor = float(self.recrossingFactor)
        k_RPMD = k_QTST * recrossingFactor
                
        fromAtomicUnits = 1e6 * ((5.2917721092e-11)**3 / 2.418884326505e-17) * constants.Na
        logging.info('Final value of rate coefficient = {0:g} cm^3/(mol*s)'.format(k_RPMD * fromAtomicUnits))
        logging.info('')

        self.saveRateCoefficient(rateFilename, k_QTST_s0, staticFactor, k_QTST, recrossingFactor, k_RPMD)

        return k_RPMD
    
    def sampleMomentum(self, Nbeads=None):
        """
        Return a pseudo-random sampling of momenta from a Boltzmann 
        distribution at the temperature of interest.
        """
        if Nbeads is None: Nbeads = self.Nbeads
        return system.sample_momentum(self.mass, self.beta, Nbeads)

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

    def initializeRandomNumberGenerator(self):
        """
        Initialize the random number generator. If a valid value is found in
        the ``randomSeed`` attribute, that value is used as the seed; otherwise
        a default value is used.
        """
        if self.randomSeed is not None and isinstance(self.randomSeed, int):
            random_init_seed(self.randomSeed)
        else:
            random_init()
