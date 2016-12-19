import numpy as np
from scipy import optimize
import sys

from astropy import modeling
from itertools import chain
from tardisnuclear.ejecta import Ejecta
from tardisnuclear.rad_trans import SimpleLateTime
from tardisnuclear.nuclear_data import DecayRadiation

from scipy import stats
from collections import OrderedDict
import pandas as pd

from astropy import units as u
import pymultinest

msun_to_cgs = u.Msun.to(u.g)
mpc_to_cm = u.Mpc.to(u.cm)


class BaseModel(modeling.Model):

    def __call__(self, *inputs, **kwargs):
        parameters = self._param_sets(raw=True)
        return self.evaluate(*chain(inputs, parameters))


class BolometricLightCurveModel(BaseModel):
    pass

class BolometricLightCurveModelIa(object):



    def __init__(self, epochs, lum_dens, lum_dens_err, ni56, ni57, co55, ti44):
        self.epochs = epochs
        self.lum_dens = lum_dens
        self.lum_dens_err = lum_dens_err
        self.ejecta = Ejecta.from_masses(Ni56=ni56 * u.Msun, Ni57=ni57 * u.Msun,
                                         Co55=co55 * u.Msun, Ti44=ti44 * u.Msun)
        self.nuclear_data = DecayRadiation(self.ejecta.get_all_children_nuc_name())
        self.rad_trans = SimpleLateTime(self.ejecta, self.nuclear_data)

    def calculate_light_curve(self, ni56, ni57, co55, ti44, fraction=1.0,
                              distance=6.4, epochs=None):

        if epochs is None:
            epochs = self.epochs
        total_mass = ni56 + ni57 + co55 + ti44

        self.ejecta.mass_g = total_mass * msun_to_cgs
        self.ejecta['Ni56'] = ni56 / total_mass
        self.ejecta['Ni57'] = ni57 / total_mass
        self.ejecta['Co55'] = co55 / total_mass
        self.ejecta['Ti44'] = ti44 / total_mass
        luminosity_density =  self.rad_trans.total_bolometric_light_curve(epochs)
        return (luminosity_density * fraction /
                (4 * np.pi * (distance * mpc_to_cm)**2))


    def calculate_individual_light_curve(self, ni56, ni57, co55, ti44, fraction=1.0,
                              distance=6.4, epochs=None):

        if epochs is None:
            epochs = self.epochs
        total_mass = ni56 + ni57 + co55 + ti44

        self.ejecta.mass_g = total_mass * msun_to_cgs
        self.ejecta['Ni56'] = ni56 / total_mass
        self.ejecta['Ni57'] = ni57 / total_mass
        self.ejecta['Co55'] = co55 / total_mass
        self.ejecta['Ti44'] = ti44 / total_mass
        luminosity_density =  self.rad_trans.bolometric_light_curve(epochs)
        return (luminosity_density * fraction /
                (4 * np.pi * (distance * mpc_to_cm)**2))


    def fitness_function(self, ni56, ni57, co55, ti44, fraction, distance):

        model_light_curve = self.calculate_light_curve(ni56, ni57, co55, ti44,
                                                 fraction, distance)
        return (model_light_curve.value - self.lum_dens)/self.lum_dens_err


    def log_likelihood(self, model_param, ndim, nparam):
        #return -5

        model_param = [model_param[i] for i in xrange(6)]
        return (-0.5 * self.fitness_function(*model_param)**2).sum()

    def simple_fit(self, ni56, ni57, co55, ti44, method='Nelder-Mead'):
        def fit_func(isotopes):
            ni57, co55, ti44 = np.abs(isotopes)
            mdl = self.evaluate(ni56, ni57, co55, ti44)
            mdl *= np.mean(self.luminosity / mdl.value)
            return ((mdl.value - self.luminosity)**2).sum()

        fit = optimize.minimize(fit_func, (ni57, co55, ti44),
                                 method=method)
        mdl = self.evaluate(ni56, *fit.x)
        norm_factor = np.mean(self.luminosity / mdl.value)
        mdl *= norm_factor

        return fit, norm_factor, mdl


    def multinest_fit(self, priors, **kwargs):


        mn_fit = pymultinest.run(self.log_likelihood, priors.prior_transform, 6,
                                 outputfiles_basename='sn11fe/fit', **kwargs)

        return mn_fit




class MultiNestResult():


    @classmethod
    def from_multinest_basename(cls, basename, parameter_names):
        """
        Reading a MultiNest result from a basename

        Parameters
        ----------

        basename: str
            basename (path + prefix) for a multinest run

        Returns
            : ~MultinestResult
        """

        posterior_data = cls.read_posterior_data(basename, parameter_names)

        return cls(posterior_data)

    @classmethod
    def from_hdf5(cls, h5_fname, key):
        """
        Reading a Multinest result from its generated HDF5 file

        Parameters
        ----------

        h5_fname: ~str
            HDF5 filename

        key: ~str
            group identifier in the store
        """

        posterior_data = pd.read_hdf(h5_fname, key)

        return cls(posterior_data)


    @staticmethod
    def read_posterior_data(basename, parameter_names):
        """
        Reading the posterior data into a pandas dataframe

        """
        posterior_data = pd.read_csv('{0}/fit.txt'.format(basename),
                           delim_whitespace=True,
                           names=['posterior', 'x'] + parameter_names)
        posterior_data.index = np.arange(len(posterior_data))
        return posterior_data

    def __init__(self, posterior_data):
        self.posterior_data = posterior_data
        self.parameter_names = [col_name for col_name in posterior_data.columns
                                if col_name not in ['x', 'posterior']]

    def calculate_sigmas(self, sigma):
        sigmas = OrderedDict()
        for parameter_name in self.parameter_names:
            posterior_data = self.posterior_data.sort(parameter_name)
            parameter_values, posterior_values = (posterior_data[parameter_name],
                                                  posterior_data['posterior'])
            posterior_cumsum = posterior_values.cumsum()

            norm_distr = stats.norm(loc=0.0, scale=1.)

            sigma_low = np.interp(norm_distr.cdf(-sigma), posterior_cumsum,
                                  parameter_values)

            sigma_high = np.interp(norm_distr.cdf(sigma), posterior_cumsum,
                                  parameter_values)


            sigmas[parameter_name] = (sigma_low, sigma_high)

        return sigmas

    @property
    def mean(self):
        if not hasattr(self, '_mean'):
            _mean = OrderedDict([(param_name,
                                  np.average(self.posterior_data[param_name],
                                             weights=
                                             self.posterior_data['posterior']))
                                 for param_name in self.parameter_names])
            self._mean = _mean

        return self._mean
