import urllib.request, urllib.error, urllib.parse
import os
import logging

logger = logging.getLogger(__name__)

import bs4
import pandas as pd
from pyne import ensdf

from astropy import units as u

#getting the data_path

from tardisnuclear.config import get_data_dir

TARDISNUCLEAR_DATA_DIR = get_data_dir()

from pyne import nucname
from abc import ABCMeta


class ENSDFReader:

    def __init__(self, fname):

        #### Still seemingly missing several bits of data from NNDC ###
        raw_decay_list = ensdf.decays(fname)[0]

        self.parent_nuc_id = raw_decay_list[0]
        self.daughter_nuc_id = raw_decay_list[1]
        self.reaction_id = raw_decay_list[2]
        self.half_life = raw_decay_list[3]  # in s
        self.half_life_uncert = raw_decay_list[4]  # in s
        self.branching_ratio = raw_decay_list[5]  # percent
        self.conversion_factor_gamma = raw_decay_list[6]
        self.conversion_factor_gamma_uncert = raw_decay_list[7]
        # for gamma intensity to photons per 100 decays of the parent

        self.conversion_factor_lepton = raw_decay_list[8]
        self.conversion_factor_lepton_uncert = raw_decay_list[9]
        # Conversion factor for electron capture / beta intensity to
        # electron captures / betas per 100 decays of the parent

        self.gamma = self._gamma_to_dataframe(raw_decay_list[11])
        self.alpha = self._alpha_to_dataframe(raw_decay_list[12])
        self.beta_minus = self._beta_minus_to_dataframe(raw_decay_list[13])
        self.beta_plus = self._beta_minus_to_dataframe(raw_decay_list[14])

    @staticmethod
    def _gamma_to_dataframe(gamma_list):
        columns = ['level_1_id', 'level_2_id', 'parent_id', 'daughter_id',
                   'energy', 'energy_uncert', 'intensity', 'intensity_uncert',
                   'electron_conversion_intensity', 'electron_conversion_intensity_uncert',
                   'total_transition_intensity', 'total_transition_intensity_uncert']
        columns += ['{0}_electron_conversion_intensity'.format(item) for item in 'klm']

        return pd.DataFrame(data=gamma_list, columns=columns)

    @staticmethod
    def _alpha_to_dataframe(alpha_list):
        if alpha_list == []:
            return None
        else:
            raise NotImplementedError('Not implemented yet')

        # parent nuclide id in state_id form
        # child nuclide id in state_id form
        # alpha energy
        # alpha intensity in percent of total alphas

    @staticmethod
    def _beta_minus_to_dataframe(beta_minus_list):
        if beta_minus_list == []:
            return None
        columns = ['parent_id', 'child_id', 'endpoint_energy', 'average_energy', 'intensity']
        return pd.DataFrame(beta_minus_list, columns=columns)

    @staticmethod
    def _beta_plus_to_dataframe(beta_plus_list):
        if beta_plus_list == []:
            return None
        columns = ['parent_id', 'child_id', 'endpoint_energy', 'average_energy', 'intensity',
                   'electron_capture_intensity']
        columns += ['{0}_electron_conversion_intensity'.format(item) for item in 'klm']
        return pd.DataFrame(beta_plus_list, columns=columns)


def _get_nuclear_database_path():
    if not os.path.exists(TARDISNUCLEAR_DATA_DIR):
        os.mkdir(TARDISNUCLEAR_DATA_DIR)
    return os.path.join(TARDISNUCLEAR_DATA_DIR, 'decay_radiation.h5')

def _sanitize_nuclear_string(nuclear_string):
    try:
        sanitized_nuclear_string = nucname.name(nuclear_string)
    except RuntimeError:
        raise ValueError('{0} not a valid isotope string'.format(
            nuclear_string))
    else:
        return sanitized_nuclear_string

def download_decay_radiation(nuclear_string):
    nuclear_string = _sanitize_nuclear_string(nuclear_string)
    base_url = ('http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?'
                'nuc={nucname}&unc=nds')
    data_url = base_url.format(nucname=nuclear_string.upper())
    logger.info('Downloading data from {0}'.format(data_url))
    nuclear_bs = bs4.BeautifulSoup(urllib.request.urlopen(data_url))

    data_list = nuclear_bs.find_all('u')
    if data_list == []:
        raise ValueError('{0} is stable and does not have decay '
                         'radiation'.format(nuclear_string))

    data_sets = []
    current_data_set = {}
    for item in data_list:
        data_name = item.get_text()
        if data_name.startswith('Dataset'):
            if current_data_set != {}:
                data_sets.append(current_data_set)
            current_data_set = {}


        if data_name in decay_radiation_parsers:
            next_table = item.find_next('table')
            data_result = decay_radiation_parsers[data_name].parse(next_table)
            current_data_set.update(data_result)

        elif data_name.startswith('Author'):
            pass
        elif data_name.startswith('Citation'):
            pass
        else:
            print("Data \"{0}\" is not recognized".format(item.get_text()))

    if current_data_set != {}:
            data_sets.append(current_data_set)
    return data_sets

def store_decay_radiation_from_ejecta(ejecta, force_update=False):
    """
    Check if all isotopes of a given ejecta are in the database and
    download if necessary.

    :param ejecta:
    :param force_update:
    :return:
    """

    for isotope in ejecta.get_all_children_nuc_name():
        print("Working on isotope", isotope)
        try:
            store_decay_radiation(isotope, force_update=force_update)
        except IOError as e:
            print(str(e))
            print("skipping")



def store_decay_radiation(nuclear_string, force_update=False):
    nuclear_string = _sanitize_nuclear_string(nuclear_string)
    fname = _get_nuclear_database_path()
    with pd.HDFStore(fname, mode='a') as ds:
        if nuclear_string in ds and not force_update:
            raise IOError('{0} is already in the database '
                          '(force_update to overwrite)'.format(nuclear_string))
        try:
            data_set_list = download_decay_radiation(nuclear_string)
        except ValueError:
            print("{0} is stable - making empty dataset".format(nuclear_string))
            ds['{0}'.format(nuclear_string)] = pd.DataFrame()
        else:
            for i, data_set in enumerate(data_set_list):
                for key, value in list(data_set.items()):
                    group_str = '{0}/data_set{1}/{2}'.format(nuclear_string, i,
                                                             key)
                    print("Writing group", group_str)
                    ds[group_str] = value
        ds.flush()
        ds.close()



def get_decay_radiation(nuclear_string, data_set_idx=0):
    nuclear_string = _sanitize_nuclear_string(nuclear_string)
    fname = _get_nuclear_database_path()

    if not os.path.exists(fname):
        if (not input('{0} not in database - download [Y/n]'.format(
                nuclear_string)).lower() == 'n'):
            store_decay_radiation(nuclear_string)
        else:
            raise ValueError('{0} not in database'.format(
                        nuclear_string))


    with pd.HDFStore(fname, mode='r') as ds:
        current_keys = [key for key in list(ds.keys())
                            if key.startswith('/{0}/data_set{1:d}'.format(
                nuclear_string, data_set_idx))]
        if len(current_keys) == 0:
            if '/{0}'.format(nuclear_string) in list(ds.keys()):
                logger.debug('{0} is stable - no decay radiation available'.format(
                    nuclear_string))
                return {}
            else:
                if (not input(
                        '{0} not in database - download [Y/n]'.format(
                            nuclear_string)).lower() == 'n'):
                    ds.close()
                    store_decay_radiation(nuclear_string)
                else:
                    raise ValueError('{0} not in database'.format(
                        nuclear_string))
        data_set = {}
        for key in current_keys:
            data_set[key.split('/')[-1]] = ds[key]
    return data_set



class BaseParser(metaclass=ABCMeta):
    @staticmethod
    def _convert_html_to_df(html_table, column_names):
        df = pd.read_html(str(html_table))[0].iloc[1:]
        df.columns = column_names
        if 'type' in column_names:
            df.type[df.type.isnull()] = ''
        return df
    @staticmethod
    def _sanititze_table(df):
        if 'energy' in df.columns:
            df.energy = df.energy.apply(lambda x: u.Quantity(
                float(x.split()[0]), u.keV).to(u.erg).value)
        if 'end_point_energy' in df.columns:
            df.end_point_energy = df.end_point_energy.apply(lambda x: u.Quantity(
                float(x.split()[0]), u.keV).to(u.erg).value)

        if 'intensity' in df.columns:
            df.intensity = df.intensity.apply(lambda x:
                                              float(x.split('%')[0])/100.)

        if 'dose' in df.columns:
            del df['dose']

        return df

    def _default_parse(self, html_table):
        df = self._convert_html_to_df(html_table, self.columns)
        df = self._sanititze_table(df)

        return {self.name : df}

    def parse(self, html_table):
        return self._default_parse(html_table)


class ElectronTableParser(BaseParser):
    html_name = 'Electrons'
    name = 'electrons'
    columns = ['type', 'energy', 'intensity', 'dose']


class BetaPlusTableParser(BaseParser):
    html_name = 'Beta+'
    name = 'beta_plus'
    columns = ['energy', 'end_point_energy', 'intensity', 'dose']


class BetaMinusTableParser(BaseParser):
    html_name = 'Beta-'
    name = 'beta_minus'
    columns = ['energy', 'end_point_energy', 'intensity', 'dose']



class XGammaRayParser(BaseParser):
    html_name = 'Gamma and X-ray radiation'
    columns = ['type', 'energy', 'intensity', 'dose']

    def parse(self, html_table):
        df = self._convert_html_to_df(html_table, self.columns)
        df = self._sanititze_table(df)

        x_ray_mask = df.type.str.startswith('XR')

        results = {}
        results['x_rays'] = df[x_ray_mask]
        results['gamma_rays'] = df[~x_ray_mask]
        return results






decay_radiation_parsers = {item.html_name:item()
                          for item in BaseParser.__subclasses__()}

