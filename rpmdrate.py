#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RPMDrate - Ring polymer molecular dynamics simulations
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

"""
This is the main execution script for RPMD, a tool for using ring polymer
molecular dynamics simulations to compute properties of chemical systems, 
including, but not limited to, chemical reaction rates. To use, simply pass the
path to an RPMD input file on the command line, e.g. ::

$ python rpmdrate.py examples/H+CH4/input.py 1000 1

You can also explicitly set the number of processors to use using the ``-p``
flag, e.g. ::

$ python rpmdrate.py examples/H+CH4/input.py 1000 1 -p 8

Only one processor is used by default.
"""

import os.path
import sys
import time
import argparse
import logging

################################################################################

def parseCommandLineArguments():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('file', metavar='FILE', type=str, nargs=1, help='a file describing the job to execute')
    parser.add_argument('T', metavar='TEMP', type=float, nargs=1, help='the temperature in K')
    parser.add_argument('Nbeads', metavar='BEADS', type=int, nargs=1, help='the number of beads')
    parser.add_argument('-p', '--processes', metavar='PROC', type=int, nargs=1, default=[1], help='the number of processors to use')

    # Options for controlling the amount of information printed to the console
    # By default a moderate level of information is printed; you can either
    # ask for less (quiet), more (verbose), or much more (debug)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-q', '--quiet', action='store_const', const=logging.WARNING, default=logging.INFO, dest='verbose', help='only print warnings and errors')
    group.add_argument('-v', '--verbose', action='store_const', const=logging.DEBUG, default=logging.INFO, dest='verbose', help='print more verbose output')
    group.add_argument('-d', '--debug', action='store_const', const=0, default=logging.INFO, dest='verbose', help='print debug information')

    return parser.parse_args()

################################################################################

def initializeLog(verbose, logFile=None):
    """
    Set up a logger for PyRate to use to print output to stdout. The
    `verbose` parameter is an integer specifying the amount of log text seen
    at the console; the levels correspond to those of the :data:`logging` module.
    """
    # Create logger
    logger = logging.getLogger()
    logger.setLevel(verbose)

    # Use custom level names for cleaner log output
    logging.addLevelName(logging.CRITICAL, 'Critical: ')
    logging.addLevelName(logging.ERROR, 'Error: ')
    logging.addLevelName(logging.WARNING, 'Warning: ')
    logging.addLevelName(logging.INFO, '')
    logging.addLevelName(logging.DEBUG, '')
    logging.addLevelName(0, '')

    # Create formatter and add to handlers
    formatter = logging.Formatter('%(levelname)s%(message)s')
    
    # Remove old handlers before adding ours
    while logger.handlers:
        logger.removeHandler(logger.handlers[0])
   
    # Create console handler; send everything to stdout rather than stderr
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(verbose)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Create file handler; always be at least verbose in the file
    if logFile:
        fh = logging.FileHandler(filename=logFile)
        fh.setLevel(min(logging.DEBUG,verbose))
        fh.setFormatter(formatter)
        logger.addHandler(fh)

################################################################################

def logHeader(level=logging.INFO):
    """
    Output a header containing identifying information about RPMDrate to the
    log.
    """
    logging.log(level, 'RPMDrate execution initiated at {0}'.format(time.asctime()))
    logging.log(level, '')

    logging.log(level, '###############################################################')
    logging.log(level, '#                                                             #')
    logging.log(level, '#                          ********                           #')
    logging.log(level, '#                          RPMDrate                           #')
    logging.log(level, '#                          ********                           #')
    logging.log(level, '#                                                             #')
    logging.log(level, '#   Version: 0.1.0 (1 May 2012)                               #')
    logging.log(level, '#                                                             #')
    logging.log(level, '#   Authors: Joshua W. Allen <jwallen@mit.edu>                #')
    logging.log(level, '#            Yury V. Suleimanov <ysuleyma@mit.edu>            #')
    logging.log(level, '#            William H. Green <whgreen@mit.edu>               #')
    logging.log(level, '#                                                             #')
    logging.log(level, '#   Website: http://github.com/GreenGroup/RPMDrate            #')
    logging.log(level, '#                                                             #')
    logging.log(level, '###############################################################')
    logging.log(level, '')

################################################################################

def logFooter(level=logging.INFO):
    """
    Output a footer to the log.
    """
    logging.log(level, '')
    logging.log(level, 'RPMDrate execution terminated at {0}'.format(time.asctime()))
    
################################################################################

if __name__ == '__main__':
    
    # Parse the command-line arguments (requires the argparse module)
    args = parseCommandLineArguments()
    
    # Determine the output directory
    outputDirectory = os.path.dirname(os.path.abspath(args.file[0]))
    
    # Initialize the logging system (both to the console and to a file in the
    # output directory)
    initializeLog(args.verbose, os.path.join(outputDirectory, 'pyrate.log'))

    # Print some information to the beginning of the log
    logHeader()
    
    # Load the input file for the job
    from rpmdrate.input import loadInputFile
    system, jobList = loadInputFile(args.file[0], args.T[0], args.Nbeads[0], args.processes[0])
    logging.info('')
    
    # Run the requested jobs
    for job, params in jobList:
        if job == 'configurations':
            dt, evolutionTime, xi_list, kforce = params
            system.generateUmbrellaConfigurations(dt, evolutionTime, xi_list, kforce)
        elif job == 'umbrella':
            dt, windows, saveTrajectories = params
            system.conductUmbrellaSampling(dt, windows, saveTrajectories)
        elif job == 'PMF':
            windows, xi_min, xi_max, bins = params
            system.computePotentialOfMeanForce(windows, xi_min, xi_max, bins)
        elif job == 'recrossing':
            dt, equilibrationTime, childTrajectories, childSamplingTime, childrenPerSampling, childEvolutionTime, xi_current, saveParentTrajectory, saveChildTrajectories = params
            system.computeRecrossingFactor(dt, equilibrationTime, childTrajectories, childSamplingTime, childrenPerSampling, childEvolutionTime, xi_current, saveParentTrajectory, saveChildTrajectories)
        elif job == 'rate':
            system.computeRateCoefficient()
    
    # Print some information to the end of the log
    logFooter()
