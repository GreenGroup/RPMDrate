#!/usr/bin/env python
# encoding: utf-8

from PES import get_potential

################################################################################

label = 'H + C2H6 -> H2 + C2H5'

reactants(
    atoms = ['C', 'H', 'H', 'H', 'C', 'H', 'H', 'H', 'H'],
    reactant1Atoms = [1,2,3,4,5,6,7,8],
    reactant2Atoms = [9],
    Rinf = (15 * 0.52918,"angstrom"),
)

transitionState(
    geometry = (
        [[-0.18458289,    1.09944478,   -0.76346873],
         [-0.09371221,    0.02584027,   -0.68919968],
         [0.60997447,    1.48587575,   -1.38455625],
         [-1.13897537,    1.35045534,   -1.20224965],
         [-0.08879895,    1.70522604,    0.60541932],
         [-0.16677838,    2.77846596,    0.51335803],
         [-0.81879013,    1.34592969,    1.31564395],
         [1.09626535,    1.34842526,    1.16925957],
         [1.89966923,    1.10653532,    1.55151009]],
        "angstrom"),
    formingBonds = [(8,9)], 
    breakingBonds = [(5,8)],
)
equivalentTransitionState(
    formingBonds=[(6,9)], 
    breakingBonds=[(5,6)],
)
equivalentTransitionState(
    formingBonds=[(7,9)], 
    breakingBonds=[(5,7)],
)

thermostat('Andersen')

#################################################################################

xi_list = numpy.arange(-0.05, 1.01, 0.01)
#xi_list = [1.0]

generateUmbrellaConfigurations(
    dt = (0.0001,"ps"),
    evolutionTime = (5,"ps"),
    xi_list = xi_list,
    kforce = 0.1 * T,
)

xi_list = numpy.arange(-0.05, 1.01, 0.01)
#xi_list = [0.22]

windows = []
for xi in xi_list:
    window = Window(xi=xi, kforce=0.1*T, trajectories=150, equilibrationTime=(20,"ps"), evolutionTime=(100,"ps"))
    windows.append(window)

conductUmbrellaSampling(
    dt = (0.0001,"ps"),
    windows = windows,
)

computePotentialOfMeanForce(windows=windows, xi_min=-0.02, xi_max=1.01, bins=5000)

computeRecrossingFactor(
    dt = (0.0001,"ps"),
    equilibrationTime = (20,"ps"),
    childTrajectories = 100000,
    childSamplingTime = (2,"ps"),
    childrenPerSampling = 100,
    childEvolutionTime = (0.05,"ps"),
)

computeRateCoefficient()
