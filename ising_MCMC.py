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

KBOLTZMANN = 1.3806488e-23      # joule/kelvin
# note 0.0019872041 kcal/(mol K)
import numpy as np

class IsingSystem:
    """ 
    the (2D) ising system
    """
    def __init__(self, size, sweeps, temperature):
        self.system = np.random.random_integers(0, 1, [size, size])
        self.system *= 2
        self.system -= 1

        self.coupling_constant = 1 # let J be 1...

        # Used to calculate delta E from neighboring interactions:
        self.neighbor_positions = np.array([[0,1], [1,0], [0,-1], [-1,0]])

        times = np.arange(sweeps)
        magnetization_timeseries = np.zeros(sweeps)
        energy_timeseries = np.zeros(sweeps)
        self.system_length = size
        self.steps_per_sweep = self.system.size
        self.beta = KBOLTZMANN * temperature
        return

    def calculate_energy(self):
        """
        The full energy function involves calculation of all pairwise
        energy interactions. We sum across rows down columns separately
        """
        energy = 0.
        for row in self.system:
            energy += np.sum(
                [self.coupling_constant * i * j for i,j in zip(row[:-1], row[1:])]
                )
        for col in np.rollaxis(self.system, -1):
            energy += np.sum(
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
                neighboring_spin = self.system[position + offset]
            except IndexError:
                # our position is some edge spin and no PBC in effect
                continue
            neighboring_spins[i] = neighboring_spin
        new_position_spin = -1 * self.system[position]
        deltaE = (np.sum(new_position_spin * neighboring_spins) -
                  np.sum(self.system[position] * neighboring_spins))
        return deltaE
            
    def calculate_magnetization(self):
        return np.sum(self.system)

    def run_simulation(self, sweeps):
        """ where steps is number of sweeps"""
        if not self.initialized:
            return
        for i in xrange(sweeps):
            for j in xrange(self.steps_per_sweep):
                selected_spin = np.random.random_integers(0, self.system_length, 2)
                deltaE = self.calculate_deltaE(selected_spin)
                if deltaE < 0:
                    self.system[selected_spin] *= -1
                    continue
                else:
                    weight = boltmann_weight(deltaE, self.beta)
                    if np.random.rand() < weight:
                        self.system[selected_spin] *= -1
            self.magnetization_timeseries[i] = self.calculate_magnetization()
            self.energy_timeseries[i] = self.calculate_energy()
        return

def boltzmann_weight(energy, beta):
    return np.exp(-1 * energy * beta)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('grid_size', help="NxN")
    parser.add_argument('temperature', type=float)
    parser.add_argument('--num-sweeps', type=int, default=2500)
    parser.add_argument('--out-pkl', required=True, help='pkl log file')
    return parser.parse_args()

def main(args):
    ising_system = IsingSystem(args.grid_size)
    ising_system.initialize_simulation(args.num_sweeps, args.temperature)
    


if __name__ == "__main__":
    main(parse_args())
