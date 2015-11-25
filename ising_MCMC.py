#!/usr/local/bin/python
"""
ising_MCMC.py

ising hamiltonian: -J sum pairs - H 

another implementation of the 2D array would be of course as 1D with
appropriate indicing

questions:
what's the coupling constant?
what's the critical temperature?
do we do temperature in beta or in K

@author - Victor Zhao

"""

import numpy as np

class IsingSystem:
    """ the (2D) ising system"""

    def __init__(self, size):
        self.system = np.random.random_integers(0, 1, [size, size])
        self.system *= 2
        self.system -= 1

        self.coupling_constant = 1 # let J be 1...

        # Used to calculate delta E from neighboring interactions:
        self.neighbor_positions = np.array([[0,1], [1,0], [0,-1], [-1,0]])

        self.initialized = False
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
        neighboring_spins = []
        for offset in self.neighbor_positions:
            try:
                neighboring_spin = self.system[position + offset]
            except IndexError:
                continue            # our position is some edge spin
            neighboring_spins.append(neighboring_spin)
        neighboring_spins = np.array(neighboring_spins)
        position_spin = -1 * self.system[position]
        return np.sum(2 * position_spin * neighboring_spins)
            
    def calculate_magnetization(self):
        return np.sum(self.system)

    def initialize_simulation(self):
        """ 
        call before running sim 
        """
        self.initialized = True
        


def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('grid_size', help="NxN")
    parser.add_argument('temperature', type=float)
    parser.add_argument('--num-sweeps', type=int, default=2500)
    return parser.parse_args()

def main(args):
    ising_system = IsingSystem(args.grid_size)

    


if __name__ == "__main__":
    main(parse_args())
