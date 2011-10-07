from astropy import units as u
from pyne.material import Material
from pyne import data

import numpy as np

from astropy import units as u

msun_to_cgs = u.Msun.to(u.g)

class Ejecta(object):
    """
    Radioactive Ejecta composition

    Parameters
    ----------

    mass: ~float
        mass in solar masses

    composition: ~dict
        a composition dictionary, e.g. {'Co56':0.5, 'Ni56':0.5}
        will be normalized to 1
    """

    @classmethod
    def from_masses(cls, **kwargs):
        """
        Initialize the ejecta from masses

        Parameters
        ----------

        **kwargs: key, value pairs
            like Co56=1e33*u.g

        """

        mass = sum(kwargs.values()).to(u.g)
        composition = {key:float(value / mass) for key, value in kwargs.items()}


        return cls(mass.value, composition)

    def __init__(self, mass, composition):
        self.mass = mass * msun_to_cgs
        self.material = Material(self._normalize_composition(composition))


    def __getitem__(self, item):
        return self.material.__getitem__(item)

    def __setitem__(self, key, value):
        self.material.__setitem__(key, value)


    @staticmethod
    def _normalize_composition(composition):
        composition_sum = np.sum(composition.values())
        normed_composition = {key:value/composition_sum
                              for key, value in composition.items()}
        return normed_composition


    def decay(self, time):
        new_material= self.material.decay(time.to(u.s).value)
        return self.__class__(self.mass, new_material)

    def __repr__(self):
        return self.material.__str__()

    @property
    def N(self):
        N = {}
        for nuc in self.material:
            mass = self.material[nuc] * self.mass
            N[nuc] =  (1 / (data.atomic_mass(nuc) * u.u) * mass).to(1)
        return N



