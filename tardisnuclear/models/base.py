import numpy as np

from astropy.modeling import FittableModel, Parameter
from astropy import units as u
from pyne import nucname

from tardisnuclear.ejecta import Ejecta, msun_to_cgs
from tardisnuclear.rad_trans import SimpleLateTime
from tardisnuclear.nuclear_data import NuclearData

mpc_to_cm = u.Mpc.to(u.cm)

class BaseBolometricLightCurve(FittableModel):
    inputs = ('epoch', )
    outputs = ('luminosity', )

    def __init__(self, **kwargs):
        super(BaseBolometricLightCurve, self).__init__(**kwargs)
        self._init_ejecta(kwargs)
        nuclear_data = NuclearData(self.ejecta.get_all_children_nuc_name())
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

    def evaluate(self, epoch, *args):
        self._update_ejecta(args)
        return self.rad_trans.total_bolometric_light_curve(epoch)


def make_bolometric_model(**kwargs):
    """
    Make a bolometric lightcurve model
    :param kwargs:
    :return:
    """
    class_dict = {}
    class_dict['__init__'] = BaseBolometricLightCurve.__init__
    init_kwargs = {}
    for isotope_name in kwargs:
        if not nucname.isnuclide(isotope_name):
            raise ValueError('{0} is not a nuclide name')
        class_dict[isotope_name.lower()] = Parameter()
        init_kwargs[isotope_name.lower()] = kwargs[isotope_name]
    BolometricLightCurve = type('BolometricLightCurve',
                                (BaseBolometricLightCurve, ), class_dict)

    return BolometricLightCurve(**init_kwargs)

class RSquared(FittableModel):
    inputs = ('luminosity', )
    outputs = ('luminosity_density', )

    #Distance in Mpc
    distance = Parameter()


    def evaluate(self, luminosity, distance):
        return luminosity / (4 * np.pi * (distance * mpc_to_cm)**2)


class BolometricChi2Likelihood(FittableModel):

    inputs = ('luminosity_density', )
    outputs = ('log_likelihood')

    def __init__(self, epochs, ):
        pass