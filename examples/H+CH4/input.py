#!/usr/bin/env python
# encoding: utf-8

from PES import get_potential, initialize_potential

################################################################################

label = 'H + CH4 -> H2 + CH3'

reactants( 
    atoms = ['H', 'C', 'H', 'H', 'H', 'H'],
    reactant1Atoms = [1,2,3,4,5],
    reactant2Atoms = [6],
    Rinf = (30 * 0.52918,"angstrom"),
)

transitionState(
    geometry = (
        [[-4.68413503,   -0.43825460,   -0.07250839],
         [-5.04748906,    0.58950601,   -0.07250840],
         [-4.68411607,    1.10337961,    0.81755453],
         [-4.58404767,    1.24489401,   -1.20768359],
         [-6.13758906,    0.58951941,   -0.07250839],
         [-4.29480419,    1.61876150,   -1.94140095]],
        "angstrom",
    ),
    formingBonds = [(4,6)], 
    breakingBonds = [(2,4)],
)

equivalentTransitionState(
    formingBonds=[(5,6)], 
    breakingBonds=[(2,5)],
)
equivalentTransitionState(
    formingBonds=[(3,6)], 
    breakingBonds=[(2,3)],
)
equivalentTransitionState(
    formingBonds=[(1,6)], 
    breakingBonds=[(2,1)],
)

thermostat('Andersen')

################################################################################

xi_list = numpy.arange(-0.05, 1.05, 0.01)
    
generateUmbrellaConfigurations(
    dt = (0.0001,"ps"),
    evolutionTime = (5,"ps"),
    xi_list = xi_list,
    kforce = 0.1 * T,
)

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
