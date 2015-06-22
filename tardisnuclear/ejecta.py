from collections import OrderedDict

from astropy import units as u
import pandas as pd
from pyne.material import Material
from pyne import data
from pyne import nucname

import numpy as np

from astropy import units as u

msun_to_cgs = u.Msun.to(u.g)
u_to_g = u.u.to(u.g)

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
    def from_yann_file(cls, fname):
        data = pd.read_table(fname,
                             names=['isotope', 'mass'], delim_whitespace=True)

        data = data.set_index('isotope')
        mass = data.mass.sum()
        data['norm_mass'] = data.mass / mass
        composition = data.norm_mass.to_dict()

        return cls(mass, composition)

    @classmethod
    def from_masses(cls, **kwargs):
        """
        Initialize the ejecta from masses

        Parameters
        ----------

        **kwargs: key, value pairs
            like Co56=1e33*u.g

        """

        mass = sum(kwargs.values()).to(u.Msun)
        composition = {key:float(value / mass) for key, value in kwargs.items()}


        return cls(mass.value, composition)

    def __init__(self, mass_msol, composition):
        self.mass_g = mass_msol * msun_to_cgs
        self.material = Material(self._normalize_composition(composition))
        self._pad_material()
        atomic_masses = self.get_masses()
        self.n_per_g = [1 / atomic_masses[item]
                            for item in self.get_all_children_nuc_name()]


    @property
    def mass(self):
        return self.mass_g * u.g



    def __getitem__(self, item):
        return self.material.__getitem__(item)

    def __setitem__(self, key, value):
        self.material.__setitem__(key, value)

    def keys(self):
        return self.material.keys()

    @property
    def isotopes(self):
        return [nucname.name(id) for id in self.keys()]

    def get_decay_const(self):
        return OrderedDict((nuc_name, data.decay_const(nuc_id))
                for nuc_id, nuc_name in zip(self.get_all_children(),
                                            self.get_all_children_nuc_name()))

    def get_half_life(self):
        return [data.half_life(nuc_id) for nuc_id in self.keys()]

    def get_masses(self):
        return {nuc_name:data.atomic_mass(nuc_id) * u_to_g
                for nuc_id, nuc_name in zip(self.keys(), self.isotopes)}

    def get_all_children(self):
        children_set = set()
        def get_child(nuc_id):
            children = data.decay_children(nuc_id)
            if len(children) != 0:
                children_set.update(children)
                for child_nuc_id in children:
                    get_child(child_nuc_id)

        for nuc_id in self.material:
            children_set.add(nuc_id)
            get_child(nuc_id)

        return sorted(children_set)

    def get_all_children_nuc_name(self):
        return [nucname.name(nuc_id) for nuc_id in self.get_all_children()]


    @staticmethod
    def _normalize_composition(composition):
        composition_sum = np.sum(composition.values())
        normed_composition = {key:value/composition_sum
                              for key, value in composition.items()}
        return normed_composition

    def _pad_material(self):
        for isotope in self.get_all_children_nuc_name():
            try:
                self.material[isotope]
            except KeyError:
                self.material[isotope] = 0.0

    def decay(self, epochs):
        """
        Decay the ejecta material

        Parameters
        ----------

        epochs: numpy or quantity array

        Returns
        -------
            : ~pd.DataFrame

        """
        epochs = u.Quantity(epochs, u.day)
        isotope_children = self.get_all_children()
        columns = self.get_all_children_nuc_name()
        decayed_fractions = pd.DataFrame(index=epochs.value, columns=columns,
                                       dtype=np.float64)
        day_to_s = 24 * 3600
        for epoch, row in decayed_fractions.iterrows():
            new_material = self.material.decay(epoch * day_to_s)
            fractions = [0.0 if key not in new_material else new_material[key]
                         for key in isotope_children]
            row.values[:] = fractions
        return decayed_fractions

    def get_decayed_numbers(self, epochs):
        epochs = u.Quantity(epochs, u.day)

        fractions = self.decay_epochs(epochs)

        return fractions * self.mass_g * self.n_per_g



    def __repr__(self):
        return self.material.__str__()

    def get_numbers(self):
        N = {}
        for nuc_id, nuc_name in zip(self.keys(), self.isotopes):
            mass = self.material[nuc_id] * self.mass_g
            N[nuc_name] =  (1 / (data.atomic_mass(nuc_id) * u_to_g) * mass)
        return N


    @property
    def N(self):
        return self.get_number()


