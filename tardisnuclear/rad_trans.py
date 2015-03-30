import numpy as np

from astropy import units as u
from tardisnuclear.io import get_decay_radiation

class BaseRadiationTransfer(object):
    def __init__(self, ejecta, nuclear_data):
        self.ejecta = ejecta
        self.nuclear_data = nuclear_data
        self.decay_const = ejecta.get_decay_const()

class SimpleLateTime(BaseRadiationTransfer):
    def __init__(self, ejecta, nuclear_data):
        super(SimpleLateTime, self).__init__(ejecta, nuclear_data)



    def calculate_leptons_energy(self, time):
        current_ejecta = self.ejecta.decay(time)
        nuclear_numbers = current_ejecta.get_numbers()

        leptons_energy = {}


        for isotope in current_ejecta.isotopes:
            if self.nuclear_data[isotope] == {}:
                continue
            l_energy_p_decay = self.nuclear_data[isotope][
                'total_lepton_energy_per_decay']
            leptons_energy[isotope] = (l_energy_p_decay *
                                       nuclear_numbers[isotope] *
                                       self.decay_const[isotope])
        return leptons_energy


    def bolometric_light_curve(self, epochs):
        energy = []

        for time in epochs:
            leptons_energy = self.calculate_leptons_energy(time)
            energy.append(np.sum(leptons_energy.values()))

        return energy * u.erg / u.s



