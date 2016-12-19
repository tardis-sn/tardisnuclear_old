import numpy as np
import pandas as pd

from astropy import units as u, constants as const
from astropy import modeling


class BaseEnergyInjection(modeling.Model):

    inputs = ()
    outputs = ('lepton', 'em')

    def __init__(self, ejecta, decay_radiation, cutoff_em_energy):
        super(BaseEnergyInjection, self).__init__()
        self.ejecta = ejecta
        self.decay_radiation = decay_radiation
        decay_constant = ejecta.get_decay_constant()
        self.decay_constant = pd.DataFrame(data=[decay_constant.values()],
                                           columns=decay_constant.keys())

        self.em_energy_per_decay = self._get_em_energy_per_decay(
            cutoff_energy=cutoff_em_energy)
        self.lepton_energy_per_decay = self._get_lepton_energy_per_decay()

    def evaluate(self):
        pass

    def _get_lepton_energy_per_decay(self):
        """
        Get the lepton energy for decay for each of the isotopes

        Returns
        =======
            : pandas.DataFrame

        """
        energies_per_decay = np.zeros(len(self.ejecta.isotopes))
        for channel in ['beta_plus', 'beta_minus', 'electrons']:
            energy_per_decay = []
            for isotope in self.ejecta.isotopes:
                decay_rad = self.decay_radiation[isotope].get(channel, None)
                if decay_rad is None:
                    energy_per_decay.append(0.0)
                else:
                    energy_per_decay.append(decay_rad.energy.sum())
            energies_per_decay += energy_per_decay

        return pd.DataFrame(data=[energies_per_decay], columns=self.ejecta.isotopes)

    def _get_em_energy_per_decay(self, cutoff_energy=np.inf):
        """
        Get the electromagnetic energy per decay for each of the isotopes

        Parameters
        ----------
        cutoff_energy : float or astropy.Quantity
            count energies up to this value into the energy contribution
            [default = +inf]

        Returns
        -------
            : pandas.DataFrame
        """
        cutoff_energy = u.Quantity(cutoff_energy, u.eV).to('erg').value
        energies_per_decay = []
        for isotope in self.ejecta.isotopes:
            xray = self.decay_radiation[isotope].get('x_rays', None)
            gammaray = self.decay_radiation[isotope].get('gamma_rays', None)

            if (xray is None) and (gammaray is None):
                energies_per_decay.append(0)
                continue
            em_table = pd.concat([xray, gammaray]).sort('energy')
            cutoff_energy_id = em_table.energy.searchsorted(cutoff_energy)[0]
            energies_per_decay.append(em_table.energy[:cutoff_energy_id].sum())

        return pd.DataFrame(data=[energies_per_decay],
                            columns=self.ejecta.isotopes)


    def evaluate(self, time):
        self.calculate_electron_energy(time)

    def calculate_electron_energy(self, time):
        return self.calculate_total_channel_energy(time, 'electrons')





class BaseRadiationTransfer(object):
    def __init__(self, ejecta, nuclear_data):
        self.ejecta = ejecta
        self.nuclear_data = nuclear_data
        self.decay_const = ejecta.get_decay_constant()
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

    def calculate_energy_output(self, epochs):
        numbers = self.ejecta.get_decayed_numbers(epochs)

        lepton_energy_output = self.calculate_total_channel_energy(numbers,
                                                                   'electrons')
        for channel in ['beta_plus', 'beta_minus']:
            lepton_energy_output += self.calculate_total_channel_energy(
                numbers, channel)

        x_ray_energy_output = self.calculate_total_channel_energy(numbers,
                                                                  'x_rays')
        gamma_ray_energy_output = self.calculate_total_channel_energy(
            numbers, 'gamma_rays')

        return (lepton_energy_output.sum(axis=1).values,
                x_ray_energy_output.sum(axis=1).values,
                gamma_ray_energy_output.sum(axis=1).values)