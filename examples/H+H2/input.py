#!/usr/bin/env python 
# encoding: utf-8 
  
from h3_bkmp2 import get_potential 
 
################################################################################
 
label = 'H + H2 -> HH + H' 
 
reactants( 
    atoms = ['H', 'H', 'H'],
    reactant1Atoms = [1,2], 
    reactant2Atoms = [3], 
    Rinf = (30 * 0.52918,"angstrom"), 
) 
 
transitionState( 
    geometry = ( 
        [[  0.000000,  0.000000, -1.7570],  
         [  0.000000,  0.000000,  0.000000],  
         [  0.000000,  0.000000,  1.7570]], 
        "bohr", 
    ), 
    formingBonds = [(2,3)],  
    breakingBonds = [(1,2)], 
) 
equivalentTransitionState( 
    formingBonds=[(1,3)],  
    breakingBonds=[(2,1)], 
) 
 
thermostat('GLE', A=('gle_A.txt','s^-1')) 
 
################################################################################# 
 
xi_list = numpy.arange(-0.05, 1.05, 0.01) 
#xi_list = [-0.05] 
 
generateUmbrellaConfigurations( 
    dt = (0.0001,"ps"), 
    evolutionTime = (5,"ps"), 
    xi_list = xi_list, 
    kforce = 0.1 * T, 
) 
 
xi_list = numpy.arange(-0.05, 1.05, 0.01) 
windows = [] 
for xi in xi_list: 
    window = Window(xi=xi, kforce=0.1*T, trajectories=200, equilibrationTime=(20,"ps"), evolutionTime=(100,"ps")) 
    windows.append(window) 
 
conductUmbrellaSampling( 
    dt = (0.0001,"ps"), 
    windows = windows, 
) 
 
computePotentialOfMeanForce(windows=windows, xi_min=-0.02, xi_max=1.02, bins=5000) 
 
computeRecrossingFactor( 
    dt = (0.0001,"ps"), 
    equilibrationTime = (20,"ps"), 
    childTrajectories = 100000, 
    childSamplingTime = (2,"ps"), 
    childrenPerSampling = 100, 
    childEvolutionTime = (0.05,"ps"), 
) 
 
computeRateCoefficient() 
 
