import numpy as np

from astropy import units as u

class BaseRadiationTransfer(object):
    def __init__(self, ejecta, nuclear_data):
        self.ejecta = ejecta
        self.nuclear_data = nuclear_data
        self.decay_const = ejecta.get_decay_const()
        self.energy_per_decay = self._get_energy_per_decay()


    def _get_energy_per_decay(self, channels=['x_rays', 'gamma_rays',
                                              'beta_plus', 'beta_minus',
                                              'electrons']):
        energy_per_decay = {}
        for channel in channels:
                energy_name = 'total_{0}_energy_per_decay'.format(channel)
                energy_per_decay[channel] = [
                    self.nuclear_data[item].get(energy_name, 0.0)
                    for item in self.ejecta.isotopes]

        return energy_per_decay



    def calculate_total_channel_energy(self, numbers, channel_name):
        """
        Calculate the total decay energy for a certain channel

        :param time:
        :param channel_energy_per_decay:
        :return:
        """

        energy_per_decay = self.energy_per_decay[channel_name]

        channel_energy = numbers * self.decay_const.values() * energy_per_decay

        return channel_energy


    def calculate_electron_energy(self, time):
        return self.calculate_total_channel_energy(time, 'electrons')


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



    def calculate_xray_energy(self, time):
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


    def total_bolometric_light_curve(self, epochs,
                               channels=['electrons', 'beta_plus',
                                         'beta_minus', 'x_rays']):



        return (self.bolometric_light_curve(epochs, channels).sum(axis=1).values
                * u.erg / u.s)


    def bolometric_light_curve(self, epochs,
                               channels=['electrons', 'beta_plus',
                                         'beta_minus', 'x_rays']):

        numbers = self.ejecta.get_decayed_numbers(epochs)

        energies = []
        for channel in channels:
            energies.append(self.calculate_total_channel_energy(numbers,
                                                                channel))

        total_energies = sum(energies)

        return total_energies

