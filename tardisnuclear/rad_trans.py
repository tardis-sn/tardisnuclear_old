import numpy as np

from astropy import units as u
from tardisnuclear.io import get_decay_radiation

class BaseRadiationTransfer(object):
    def __init__(self, ejecta, nuclear_data):
        self.ejecta = ejecta
        self.decay_const = ejecta.get_decay_const()
        self.decay_radiation_data = self._get_decay_radiation_data(ejecta)

    @staticmethod
    def _get_decay_radiation_data(ejecta):
        decay_radiation = {}
        for nuclear_name in ejecta.get_all_children_nuc_name():
            decay_radiation[nuclear_name] = get_decay_radiation(nuclear_name)
        return decay_radiation


class SimpleLateTime(BaseRadiationTransfer):
    def __init__(self, ejecta, nuclear_data):
        super(SimpleLateTime, self).__init__(ejecta)
        current_ejecta = self.ejecta.decay(1 * u.s)
        self.total_electron_energy_per_decay = {}
        self.total_positron_energy_per_decay = {}

        for nuclear_name in current_ejecta.isotopes:
            nuclear_rad_data = self.decay_radiation_data[nuclear_name]
            if 'beta_plus' in nuclear_rad_data:
                current_data = nuclear_rad_data['beta_plus']
                self.total_positron_energy_per_decay[nuclear_name] = (
                    current_data.energy * current_data.intensity).sum()

            if 'beta_minus' in nuclear_rad_data:
                current_data = nuclear_rad_data['beta_minus']
                self.total_electron_energy_per_decay[nuclear_name] = (
                    current_data.energy * current_data.intensity).sum()


            if 'electrons' in nuclear_rad_data:
                current_data = nuclear_rad_data['electrons']
                self.total_electron_energy_per_decay[nuclear_name] = (
                    current_data.energy * current_data.intensity).sum()




    def calculate_decay_energy(self, time):
        current_ejecta = self.ejecta.decay(time)
        nuclear_numbers = current_ejecta.get_numbers()

        electron_energy = {}
        beta_plus_energy = {}
        for nuclear_name in current_ejecta.isotopes:
            nuclear_number = nuclear_numbers[nuclear_name]

            if nuclear_name in self.total_electron_energy_per_decay:
                electron_energy[nuclear_name] = (
                    nuclear_number *
                    self.total_electron_energy_per_decay[nuclear_name]
                    * self.decay_const[nuclear_name])

            if nuclear_name in self.total_positron_energy_per_decay:
                beta_plus_energy[nuclear_name] = (
                    nuclear_number *
                    self.total_positron_energy_per_decay[nuclear_name] *
                    self.decay_const[nuclear_name])

        return electron_energy, beta_plus_energy


    def bolometric_light_curve(self, epochs):
        energy = []

        for time in epochs:
            e_energy, p_energy = self.calculate_decay_energy(time)
            energy.append(np.sum(e_energy.values() + p_energy.values()))

        return energy * u.erg / u.s



