import numpy
from amber import templates
import copy
import logging
import force_scaling


class DistRestraint:
    """ For now, this is a harmonic constraint over a squared distance D = d^2
     where D = sum_{i,j} d^2_ij over all contacts."""

    def __init__(self,contacts,kspring,force_scaler=None):
        """Initialize the DistRestraint object"""
        self.contacts = contacts        # a list of duples (start of chain is 0)
        self.kspring = kspring          # spring constant for restraint
                        # (J/[lattice space]^2)
        if force_scaler is None:
            self._force_scaler = force_scaling.NoScaling()
        else:
            self._force_scaler = force_scaler

    def evaluate_energy(self,chain,alphas):
        """Return the sum of squared-distances over the selected contacts. And transform into energies
        here the sum of distances times the force contstant is already an energy"""
        D = 0.0
        for i in range(0,len(self.contacts)):
            c = self.contacts[i][0]
            d = self.contacts[i][1]
            D = D + (chain.coords[c][0]-chain.coords[d][0])*(chain.coords[c][0]-chain.coords[d][0])
            D = D + (chain.coords[c][1]-chain.coords[d][1])*(chain.coords[c][1]-chain.coords[d][1])

        D = D * self.kspring
        scale_factors = [ self._force_scaler.scale_for_alpha(alpha) for alpha in alphas ]
        return numpy.array([D*s for s in scale_factors])

    def assign(self, alpha, weight):
        """Return the k and contacts to run the MC steps"""
        TODO
        scaled_k = self.kspring * self._force_scaler.scale_for_alpha(alpha) * weight
        contacts = se

        return scaled_k, contacts

