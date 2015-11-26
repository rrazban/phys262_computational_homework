#!/usr/local/bin/python
"""
ising_MCMC.py

ising hamiltonian: -J sum pairs - H sum spins
but i think H = 0
 
questions:
-what's the coupling constant?
-what's the critical temperature? (i think lecture notes describe it in
terms of J and zeta)
-do we do temperature in beta or in K

@author - Victor Zhao

"""

# KBOLTZMANN = 1.3806488e-23      # joule/kelvin
# note 0.0019872041 kcal/(mol K)
import numpy as np
import cPickle as pkl
import sys

class IsingSystem:
    """ 
    the (2D) ising system
    """
    def __init__(self, size, sweeps, temperature, external_field=0.,
                 coupling_constant=1.):
        self.system = np.random.random_integers(0, 1, [size, size])
        self.system *= 2
        self.system -= 1

        self.coupling_constant = coupling_constant # let J be 1
        self.external_field = external_field

        # Used to calculate delta E from neighboring interactions:
        self.neighbor_positions = np.array([[0,1], [1,0], [0,-1], [-1,0]])

        self.times = np.arange(sweeps+1)
        self.magnetization_timeseries = np.zeros(sweeps+1)
        self.energy_timeseries = np.zeros(sweeps+1)
        self.success_timeseries = np.zeros(sweeps+1, dtype=int) # track number of accepted moves per sweep
        self.system_length = size
        self.steps_per_sweep = self.system.size
        # self.beta = 1/(KBOLTZMANN * temperature)
        self.beta = 1/temperature
        return

    def calculate_energy(self):
        """
        The full energy function involves calculation of all pairwise
        energy interactions. We sum across rows down columns separately
        """
        energy = 0.
        for row in self.system:
            energy += -1 * np.sum(
                [self.coupling_constant * i * j for i,j in zip(row[:-1], row[1:])]
                )
            energy += -1 * np.sum(self.external_field * row)
        for col in np.rollaxis(self.system, -1):
            energy += -1 * np.sum(
                [self.coupling_constant * i * j for i,j in zip(col[:-1], col[1:])]
                )
        return energy
    
    def calculate_deltaE(self, position):
        """
        Position had better be a 2-tuple or list of size 2.
        """
        neighboring_spins = np.zeros(4)
        for i, offset in enumerate(self.neighbor_positions):
            try:
                neighboring_spin = self.system[tuple(position + offset)]
            except IndexError:
                # our position is some edge spin and no PBC in effect
                continue
            neighboring_spins[i] = neighboring_spin
        current_E = -1 * self.coupling_constant * np.sum(
            self.system[tuple(position)] * neighboring_spins)
        current_E += -1 * self.external_field * self.system[tuple(position)]
        new_E = -1 * current_E
        deltaE = new_E - current_E
        return deltaE
            
    def calculate_magnetization(self):
        return np.sum(self.system)

    def run_simulation(self, verbose=False):
        """ where steps is number of sweeps"""
        self.magnetization_timeseries[0] = self.calculate_magnetization()
        self.energy_timeseries[0] = self.calculate_energy()
        self.success_timeseries[0] = 0
        if verbose:
            print "{0:>10} {1:>15} {2:>15} {3:>15}".format(
                "Time", "Magnetization", "Energy", "NumSuccesses")
            print "{0:10d} {1:15.3f} {2:15.3f} {3:15d}".format(
                0, self.magnetization_timeseries[0],
                self.energy_timeseries[0], self.success_timeseries[0])
        for time in self.times[1:]:
            num_successes = 0
            for j in xrange(self.steps_per_sweep):
                selected_spin = np.random.random_integers(0, self.system_length - 1, 2)
                deltaE = self.calculate_deltaE(selected_spin)
                if deltaE < 0:
                    self.system[tuple(selected_spin)] *= -1
                    num_successes += 1
                    # print selected_spin, np.sum(self.system)
                    # print self.system
                    continue
                else:
                    weight = boltzmann_weight(deltaE, self.beta)
                    if np.random.rand() < weight:
                        self.system[tuple(selected_spin)] *= -1
                        num_successes += 1
                        # print selected_spin, np.sum(self.system)
                        # print self.system
            self.magnetization_timeseries[time] = self.calculate_magnetization()
            self.energy_timeseries[time] = self.calculate_energy()
            self.success_timeseries[time] = num_successes
            if verbose:
                print "{0:10d} {1:15.3f} {2:15.3f} {3:15d}".format(
                    time, self.magnetization_timeseries[time],
                    self.energy_timeseries[time], self.success_timeseries[time])
        return

def boltzmann_weight(energy, beta):
    return np.exp(-1 * energy * beta)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('grid_size', type=int, help="Provide N, grid will be NxN")
    parser.add_argument('temperature', type=float, help="units of kT")
    parser.add_argument('--num-sweeps', type=int, default=2500)
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--out-pkl', required=True, help='pkl log file')
    parser.add_argument('--H', type=float, default=0., help="external magnetic field")
    return parser.parse_args()

def main(args):
    print "# ising_MCMC.py"
    print "# Initializing 2D ising system of size %dx%d" % (args.grid_size, args.grid_size)
    print "# Will run at %.3f kT for %d sweeps" % (args.temperature, args.num_sweeps)

    ising_system = IsingSystem(args.grid_size, args.num_sweeps, args.temperature,
                               external_field=args.H)
    ising_system.run_simulation(args.verbose)

    print "# Simulation completed"
    with open(args.out_pkl, 'w') as f:
        pkl.dump(dict(
            times=ising_system.times,
            magnetizations=ising_system.magnetization_timeseries,
            energies=ising_system.energy_timeseries), f)

if __name__ == "__main__":
    main(parse_args())
