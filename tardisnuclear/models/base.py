import numpy as np

from astropy.modeling import FittableModel, Parameter
from astropy import units as u
from pyne import nucname
import pandas as pd

from tardisnuclear.ejecta import Ejecta, msun_to_cgs

from tardisnuclear.nuclear_data import DecayRadiation

mpc_to_cm = u.Mpc.to(u.cm)

class BaseEnergyInjection(FittableModel):
    inputs = ('time', )
    outputs = ('time', 'total')

    standard_broadcasting = False

    def __init__(self, cutoff_em_energy, **kwargs):
        super(BaseEnergyInjection, self).__init__(**kwargs)

        self._init_ejecta(kwargs)
        decay_constant = self.ejecta.get_decay_constant()

        self.decay_constant = pd.DataFrame(data=[decay_constant.values()],
                                           columns=decay_constant.keys())
        self.decay_radiation = DecayRadiation(
            self.ejecta.get_all_children_nuc_name())

        self.em_energy_per_decay = self._get_em_energy_per_decay(
            cutoff_energy=cutoff_em_energy)
        self.lepton_energy_per_decay = self._get_lepton_energy_per_decay()


    def _init_ejecta(self, isotope_dict):
        titled_isotope_dict = {name.title() : value * u.Msun
                               for name, value in isotope_dict.items()}
        self.ejecta = Ejecta.from_masses(**titled_isotope_dict)

    def _update_ejecta(self, isotope_masses):
        assert len(isotope_masses) == len(self.param_names)
        total_mass = np.sum(isotope_masses)
        self.ejecta.mass_g = total_mass * msun_to_cgs
        for isotope_name, isotope_mass in zip(self.param_names, isotope_masses):
            self.ejecta[isotope_name.title()] = isotope_mass / total_mass

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
        data = [
            self.decay_radiation[isotope].get(
                'total_lepton_energy_per_decay', 0.0)
            for isotope in self.ejecta.isotopes]
        return pd.DataFrame(data=[data], columns=self.ejecta.isotopes)

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
            em_table = pd.concat([xray, gammaray]).sort_values('energy')
            cutoff_energy_id = em_table.energy.searchsorted(cutoff_energy)
            energies_per_decay.append(
                (em_table.energy[:cutoff_energy_id] *
                 em_table.intensity[:cutoff_energy_id]).sum())

        return pd.DataFrame(data=[energies_per_decay],
                            columns=self.ejecta.isotopes)


    def calculate_injected_energy_per_s(self, time):
        energy_per_s = ((
            (self.em_energy_per_decay + self.lepton_energy_per_decay) *
            self.decay_constant).values * self.ejecta.get_decayed_numbers(time))
        return energy_per_s

    def evaluate(self, time, *args):
        self._update_ejecta(args)
        return (time, self.calculate_injected_energy_per_s(time).sum(axis=1).values)

def make_energy_injection_model(cutoff_em_energy=20*u.keV, **kwargs):
    """
    Make a bolometric lightcurve model
    :param kwargs:
    :return:
    """
    class_dict = {}
    class_dict['__init__'] = BaseEnergyInjection.__init__
    #class_dict['evaluate'] = BaseBolometricLightCurve.evaluate_specific

    init_kwargs = {}
    for isotope_name in kwargs:
        if not nucname.isnuclide(isotope_name):
            raise ValueError('{0} is not a nuclide name')
        class_dict[isotope_name.lower()] = Parameter()
        init_kwargs[isotope_name.lower()] = kwargs[isotope_name]

    EnergyInjection = type('EnergyInjection',
                                (BaseEnergyInjection,), class_dict)

    return EnergyInjection(cutoff_em_energy, **init_kwargs)

class RSquared(FittableModel):
    inputs = ('epoch', 'luminosity', )
    outputs = ('epoch', 'luminosity_density',)

    #Distance in Mpc
    distance = Parameter()


    def evaluate(self, epoch, luminosity, distance):
        luminosity_density = (luminosity /
                              (4 * np.pi * (distance * mpc_to_cm)**2))
        return epoch, luminosity_density


class BolometricChi2Likelihood(FittableModel):

    inputs = ('luminosity_density', )
    outputs = ('log_likelihood',)

    def __init__(self, epochs, ):
        pass

class SEDModel(FittableModel):
    """


    Parameters
    ----------
    spectrum : specutils.Spectrum1D
    filters : wsynphot.FilterSet
    distance : float
        distance to the supernova in Mpc
    fraction : float
        fraction of luminosity in the spectrum
    """

    inputs = ('luminosity')
    outputs = ('flux')

    distance = Parameter()
    fraction = Parameter(default=1.0)


    def __init__(self, spectrum, distance, fraction=1.0):
        assert False # unusable for now
        super(SEDModel, self).__init__(distance=distance, fraction=fraction)
        assert spectrum.unit == u.erg / u.s / u.angstrom
        luminosity_density = spectrum.data * spectrum.unit
        norm_factor = np.trapz(luminosity_density, spectrum.wavelength)
        self.normed_luminosity_density = (luminosity_density / norm_factor)

    def evaluate(self, luminosity, distance, fraction):
        luminosity_density = (
            self.normed_luminosity_density * luminosity * fraction)
        return luminosity_density / 4 * np.pi * (distance * mpc_to_cm) ** 2
