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
###############################################################################

"""
This module contains functionality for parsing RPMD input files.
"""

import os.path
import sys
import numpy
import logging

from .surface import *
from .thermostat import *
from .main import *

################################################################################

class InputError(Exception):
    """
    An exception raised when an invalid input file is read. Pass a string
    describing the exceptional behavior.
    """
    pass

################################################################################

reactants = None
transitionState = None
equivalentTransitionStates = []
thermostat = None
jobList = []

def setReactants(atoms, reactant1Atoms, reactant2Atoms, Rinf):
    global reactants
    reactants = Reactants(atoms, reactant1Atoms, reactant2Atoms, Rinf)

def setTransitionState(geometry, formingBonds, breakingBonds):
    global transitionState
    transitionState = TransitionState(geometry=geometry, formingBonds=formingBonds, breakingBonds=breakingBonds)

def addEquivalentTransitionState(formingBonds, breakingBonds):
    global equivalentTransitionStates
    equivalentTransitionStates.append([formingBonds, breakingBonds])

def setThermostat(type, **kwargs):
    global thermostat
    thermostat = [type, kwargs]

def generateUmbrellaConfigurations(dt, evolutionTime, xi_list, kforce):
    global jobList
    jobList.append(['configurations', (dt, evolutionTime, xi_list, kforce)])

def conductUmbrellaSampling(dt, windows, saveTrajectories=False):
    global jobList
    jobList.append(['umbrella', (dt, windows, saveTrajectories)])

def computePotentialOfMeanForce(windows, xi_min, xi_max, bins):
    global jobList
    jobList.append(['PMF', (windows, xi_min, xi_max, bins)])

def computeRecrossingFactor(dt, equilibrationTime, childTrajectories, childSamplingTime, childrenPerSampling, childEvolutionTime, xi_current=None, saveParentTrajectory=False, saveChildTrajectories=False):
    global jobList
    jobList.append(['recrossing', (dt, equilibrationTime, childTrajectories, childSamplingTime, childrenPerSampling, childEvolutionTime, xi_current, saveParentTrajectory, saveChildTrajectories)])

def computeRateCoefficient():
    global jobList
    jobList.append(['rate', tuple()])

getPotential = None
def potential(q):
    global getPotential
    return getPotential(q)

################################################################################

def loadInputFile(path, T, Nbeads, processes=1):
    """
    Load the RPMD input file located at `path`.
    """
    global reactants, transitionState, equivalentTransitionStates, thermostat, jobList, getPotential
    
    logging.info('Reading input file {0!r}...'.format(path))
    
    reactants = None
    transitionState = None
    equivalentTransitionStates = []
    thermostat = None
    jobList = []
    
    errorList = []
    
    global_context = {}
    local_context = {
        'reactants': setReactants,
        'transitionState': setTransitionState,
        'equivalentTransitionState': addEquivalentTransitionState,
        'thermostat': setThermostat,
        'generateUmbrellaConfigurations': generateUmbrellaConfigurations,
        'conductUmbrellaSampling': conductUmbrellaSampling,
        'computePotentialOfMeanForce': computePotentialOfMeanForce,
        'computeRecrossingFactor': computeRecrossingFactor,
        'computeRateCoefficient': computeRateCoefficient,
        'numpy': numpy,
        'T': T,
        'Nbeads': Nbeads,
        'Window': Window,
    }
    
    sys.path.append(os.path.dirname(path))
    f = open(path, 'r')
    try:
        exec f in global_context, local_context
    except (NameError, TypeError, SyntaxError), e:
        logging.error('The input file {0!r} was invalid:'.format(path))
        raise
    f.close()
    sys.path.pop()

    label = local_context.get('label', '')
    getPotential = local_context.get('get_potential', None)
    initializePotential = local_context.get('initialize_potential', None)
    
    if not getPotential:
        errorList.append('No potential energy surface supplied; you must specify a PES via the function get_potential().')
    
    if thermostat is None:
        errorList.append('No thermostat supplied; please provide a thermostat() block.')
    else:
        thermostatType, thermostatArgs = thermostat
        
        if thermostatType.lower() == 'andersen':
            samplingTime = thermostatArgs.get('samplingTime', None)
            thermostat = AndersenThermostat(samplingTime=samplingTime)
        elif thermostatType.lower() == 'gle':
            A = thermostatArgs.get('A', None)
            C = thermostatArgs.get('C', None)
            if isinstance(A,(list,tuple)) and len(A) == 2: 
                A = (os.path.join(os.path.dirname(path), A[0]), A[1])
                if not os.path.exists(A[0]):
                    errorList.append('Invalid path to A matrix for GLE thermostat.')
            if isinstance(C,(list,tuple)) and len(C) == 2: 
                C = (os.path.join(os.path.dirname(path), C[0]), C[1])
                if not os.path.exists(C[0]):
                    errorList.append('Invalid path to C matrix for GLE thermostat.')
            thermostat = GLEThermostat(A=A, C=C)
            if A is None:
                errorList.append('To use the GLE thermostat, you must specify an A matrix.')
        else:
            errorList.append('Invalid thermostat {0!r}; valid thermostats are Andersen and GLE.'.format(thermostatType))
    
    # Sanity checking of the specified breaking and forming bond pairs in the transition state
    Nbonds = transitionState.formingBonds.shape[0]
    if transitionState.breakingBonds.shape[0] != Nbonds:
        errorList.append('Different number of forming and breaking bonds encountered in transition state; these must be paired.')
    for formingBonds, breakingBonds in equivalentTransitionStates:
        if len(formingBonds) != Nbonds or len(breakingBonds) != Nbonds:
            errorList.append('Different number of forming and breaking bonds encountered in an equivalent transition state; these must correspond to those in the first transition state.')
    for n in range(Nbonds):
        atom1Found = transitionState.formingBonds[n,0] in transitionState.breakingBonds[n,:]
        atom2Found = transitionState.formingBonds[n,1] in transitionState.breakingBonds[n,:]
        if atom1Found and atom2Found:
            errorList.append('In main transition state, forming bond ({0:d},{1:d}) is the same as a breaking bond ({2:d},{3:d}).'.format(
                transitionState.formingBonds[n,0],
                transitionState.formingBonds[n,1],
                transitionState.breakingBonds[n,0],
                transitionState.breakingBonds[n,1],
            ))
        elif not atom1Found and not atom2Found:
            errorList.append('In main transition state, forming bond ({0:d},{1:d}) and breaking bond ({2:d},{3:d}) do not share an atom in common.'.format(
                transitionState.formingBonds[n,0],
                transitionState.formingBonds[n,1],
                transitionState.breakingBonds[n,0],
                transitionState.breakingBonds[n,1],
            ))
    for formingBonds, breakingBonds in equivalentTransitionStates:
        for n in range(Nbonds):
            atom1Found = formingBonds[n][0] in breakingBonds[n]
            atom2Found = formingBonds[n][1] in breakingBonds[n]
            if atom1Found and atom2Found:
                errorList.append('In an equivalent transition state, forming bond ({0:d},{1:d}) is the same as a breaking bond ({2:d},{3:d}).'.format(
                    formingBonds[n][0],
                    formingBonds[n][1],
                    breakingBonds[n][0],
                    breakingBonds[n][1],
                ))
            elif not atom1Found and not atom2Found:
                errorList.append('In an equivalent transition state, forming bond ({0:d},{1:d}) and breaking bond ({2:d},{3:d}) do not share an atom in common.'.format(
                    formingBonds[n][0],
                    formingBonds[n][1],
                    breakingBonds[n][0],
                    breakingBonds[n][1],
                ))
        
        
    if errorList:
        raise InputError('The input file {0!r} was invalid:\n- {1}'.format(path, '\n- '.join(errorList)))
    
    system = RPMD(
        label = label,
        T = T,
        Nbeads = Nbeads,
        reactants = reactants, 
        transitionState = transitionState, 
        potential = potential,
        thermostat = thermostat,
        processes = processes,
        outputDirectory = os.path.dirname(path),
    )
    for formingBonds, breakingBonds in equivalentTransitionStates:
        system.addEquivalentTransitionState(formingBonds, breakingBonds)
    
    if initializePotential:
        initializePotential()

    return system, jobList
