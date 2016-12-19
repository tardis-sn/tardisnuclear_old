import itertools
import numpy as np

from astropy.modeling import FittableModel, Parameter
from astropy import units as u
from pyne import nucname

from tardisnuclear.ejecta import Ejecta, msun_to_cgs
from tardisnuclear.rad_trans import SimpleLateTime
from tardisnuclear.nuclear_data import DecayRadiation

mpc_to_cm = u.Mpc.to(u.cm)

class BaseEnergyInjection(FittableModel):
    inputs = ('epoch', )
    outputs = ('epoch', 'lepton', 'x-ray', 'gamma-ray')

    standard_broadcasting = False

    def __init__(self, **kwargs):
        super(BaseEnergyInjection, self).__init__(**kwargs)
        self._init_ejecta(kwargs)
        nuclear_data = DecayRadiation(self.ejecta.get_all_children_nuc_name())
        self.rad_trans = SimpleLateTime(self.ejecta, nuclear_data)

    def _init_ejecta(self, isotope_dict):
        titled_isotope_dict = {name.title():value * u.Msun
                               for name, value in isotope_dict.items()}
        self.ejecta = Ejecta.from_masses(**titled_isotope_dict)

    def _update_ejecta(self, isotope_masses):
        assert len(isotope_masses) == len(self.param_names)
        total_mass = np.sum(isotope_masses)
        self.ejecta.mass_g = total_mass * msun_to_cgs
        for isotope_name, isotope_mass in zip(self.param_names, isotope_masses):
            self.ejecta[isotope_name.title()] = isotope_mass / total_mass

    def evaluate_simple_combined(self, epoch, *args):
        self._update_ejecta(args)
        return (self.rad_trans.bolometric_light_curve(
            epoch, channels=self.channels), epoch)

    def evaluate(self, epoch, *args):
        self._update_ejecta(args)
        return tuple(
            itertools.chain((epoch,),
                            self.rad_trans.calculate_energy_output(epoch)))

def make_energy_injection_model(**kwargs):
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

    return EnergyInjection(**init_kwargs)

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
