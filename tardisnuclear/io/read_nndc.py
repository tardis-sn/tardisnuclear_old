import urllib2
import os

import pandas as pd
import bs4

from astropy import units as u

import pandas as pd

import tardisnuclear

data_path = os.path.join(tardisnuclear.__path__[0], 'data')





def download_decay_radiation(nuclear_string):
    nuclear_string = nuclear_string.title()
    base_url = ('http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?'
                'nuc={nucname}&unc=nds')

    nuclear_bs = bs4.BeautifulSoup(urllib2.urlopen(
        base_url.format(nucname=nuclear_string.upper())))

    results = {}

    data_list = nuclear_bs.find_all('u')
    if data_list == []:
        raise ValueError('{0} is stable and does not have decay '
                         'radiation'.format(nuclear_string))

    data_sets = []
    current_data_set = {}
    for item in data_list:
        if item.get_text().startswith('Dataset'):
            if current_data_set != {}:
                data_sets.append(current_data_set)
            current_data_set = {}

        if item.get_text() in decay_radiation_parsers:
            next_table = item.find_next('table')
            name, data = decay_radiation_parsers[item.get_text()](
                next_table)
            current_data_set[name] = data
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

    for isotope in ejecta.isotopes:
        try:
            store_decay_radiation(isotope, force_update=force_update)
        except IOError as e:
            print str(e)
            print "skipping"



def store_decay_radiation(nuclear_string, force_update=False):
    nuclear_string = nuclear_string.title()
    fname = os.path.join(data_path, 'decay_radiation.h5')
    with pd.HDFStore(fname, mode='a') as ds:
        if nuclear_string in ds and not force_update:
            raise IOError('{0} is already in the database '
                          '(force_update to overwrite)'.format(nuclear_string))
        try:
            data_set = download_decay_radiation(nuclear_string)
        except ValueError:
            print "{0} is stable - making empty dataset".format(nuclear_string)
            ds['{0}'.format(nuclear_string)] = pd.DataFrame()
        else:
            for i, data_set in enumerate(data_set):
                for key, value in data_set.items():
                    ds['{0}/data_set{1}/{2}'.format(
                        nuclear_string, i, key)] = value


def get_decay_radiation(nuclear_string, data_set_idx=0):
    nuclear_string = nuclear_string.title()
    fname = os.path.join(data_path, 'decay_radiation.h5')
    with pd.HDFStore(fname, mode='r') as ds:
        current_keys = [key for key in ds.keys()
                            if key.startswith('/{0}/data_set{1:d}'.format(
                nuclear_string, data_set_idx))]
        if len(current_keys) == 0:

            if '/{0}'.format(nuclear_string) in ds.keys():
                print '{0} is stable - no decay radiation available'.format(
                    nuclear_string)
                return {}
            else:
                raise ValueError('{0} not in database'.format(nuclear_string))
        data_set = {}
        for key in current_keys:
            data_set[key.split('/')[-1]] = ds[key]
    return data_set




def parse_electrons_table(electrons_table):
    columns = ['type', 'energy', 'intensity', 'dose']
    df = pd.read_html(unicode(electrons_table))[0].iloc[1:]
    df.columns = columns

    df.energy = df.energy.apply(lambda x: u.Quantity(float(x.split()[0]),
                                                     u.keV).to(u.erg).value)
    df.intensity = df.intensity.apply(lambda x: float(x.split('%')[0])/100.)
    del df['dose']

    return 'electrons', df

def parse_x_gamma_table(x_gamma_table):
    columns = ['type', 'energy', 'intensity', 'dose']

    df = pd.read_html(unicode(x_gamma_table))[0].iloc[1:]
    df.columns = columns

    df.energy = df.energy.apply(lambda x: u.Quantity(float(x.split()[0]),
                                                     u.keV).to(u.erg).value)
    df.intensity = df.intensity.apply(lambda x: float(x.split('%')[0])/100.)
    del df['dose']
    return 'x_gamma_rays', df

def parse_beta_plus_table(beta_plus_table):
    columns = ['energy', 'end_point_energy', 'intensity', 'dose']

    df = pd.read_html(unicode(beta_plus_table))[0].iloc[1:]
    df.columns = columns

    df.energy = df.energy.apply(lambda x: u.Quantity(float(x.split()[0]),
                                                     u.keV).to(u.erg).value)

    df.end_point_energy = df.end_point_energy.apply(lambda x: u.Quantity(
        float(x.split()[0]), u.keV).to(u.erg).value)
    df.intensity = df.intensity.apply(lambda x: float(x.split('%')[0])/100.)
    del df['dose']

    return 'beta_plus', df


decay_radiation_parsers = {'Electrons': parse_electrons_table,
                           'Beta+': parse_beta_plus_table,
                           'Gamma and X-ray radiation': parse_x_gamma_table}