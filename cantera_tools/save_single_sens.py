from tools import CanteraSimulation
import cantera as ct
import numpy as np
from copy import deepcopy
from cPickle import dump
import sys


def main():
    simulation = CanteraSimulation(mechfile='3000_SO.cti',
                                   intial_mol_fracs={'PDD(1)': 0.25, 'Toluene(2)': 0.25, 'C11(36)': 0.5},
                                   temp_0=673,
                                   pressure_0=101325 * 300,
                                   )

    cantera_reactor = ct.IdealGasConstPressureReactor(contents=simulation.solution, energy='off')
    reactor_net = ct.ReactorNet([cantera_reactor])
    reactor_net.atol = simulation.atol
    reactor_net.rtol = simulation.rtol

    m = int(sys.argv[-1])
    species_name = simulation.solution.species_name(m)
    print 'Calculating Sensitivity of All Species to Thermochemistry of {0}'.format(species_name)
    cantera_reactor.add_sensitivity_species_enthalpy(m)

    tf = 3600 * 6
    times = np.linspace(0, tf, 101)[1:]

    sensitivity = np.zeros((simulation.solution.n_species + 2, len(times)))
    for i, t in enumerate(times):
        reactor_net.advance(t)
        sensitivity[:, i] = deepcopy(reactor_net.sensitivities()[:, 0])

    print 'Saving results to file {0}.pkl'.format(m)

    with open('{0}.pkl'.format(m), 'wb') as f:
        dump(sensitivity, f)

    print 'Job Complete for {0}'.format(species_name)


if __name__ == '__main__':
    main()
