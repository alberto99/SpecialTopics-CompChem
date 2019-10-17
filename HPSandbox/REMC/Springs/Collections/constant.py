import numpy


class ConstantCollection(object):
    def __init__(self):
        self.restraints = []

    def add_restraint(self, rest):
        self.restraints.append(rest)

    def evaluate(self, structure, temperature, current_index, alphas):
        energies = []
        for index, restraint in enumerate(self.restraints):
            energy = restraint.evaluate_energy(structure, temperature, alphas)
            energies.append(energy)
        energies = numpy.array(energies)
        return numpy.sum(energies, axis=0)

    def to_string(self, alpha, scale_factor):
        restraint_string = ''
        for index, restraint in enumerate(self.restraints):
            restraint_string += restraint.to_string(alpha, scale_factor)
        return restraint_string

    def write_weights(self, output_file):
        pass

