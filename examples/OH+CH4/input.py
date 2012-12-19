#!/usr/bin/env python
# encoding: utf-8

from PES import get_potential, initialize_potential

################################################################################

label = 'OH + CH4 -> H2O + CH3'

reactants(
    atoms = ['H', 'C', 'H', 'H', 'H', 'O', 'H'],
    reactant1Atoms = [1,2,3,4,5],
    reactant2Atoms = [6,7],
    Rinf = (15 * 0.52918,"angstrom"),
)

transitionState(
    geometry = (
        [[-4.62878267,    1.25606861,    0.95459788],
         [-4.85261637,    2.15380812,    0.37457524],
         [-4.27740626,    2.99438311,    0.76831501],
         [-4.53003946,    1.95346377,   -0.88649386],
         [-5.91912714,    2.37708958,    0.44643735],
         [-4.21407574,    1.75722671,   -2.12170961],
         [-3.93964920,    2.62885961,   -2.44704966]],
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
