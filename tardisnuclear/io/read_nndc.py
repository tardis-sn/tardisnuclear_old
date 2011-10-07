import urllib2

import pandas as pd
import bs4



def read_decay_radition(nuclear_string):
    base_url = ('http://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?'
                'nuc={nucname}&unc=nds')

    nuclear_bs = bs4.BeautifulSoup(urllib2.urlopen(
        base_url.format(nucname=nuclear_string.upper())))

    results = {}
    for item in nuclear_bs.find_all('u'):
        if item.get_text() in decay_radiation_parsers:
            next_table = item.find_next('table')
            name, data = decay_radiation_parsers[item.get_text()](
                next_table)
            results[name] = data
    return results



def parse_electrons_table(electrons_table):
    columns = ['type', 'energy', 'intensity', 'dose']
    df = pd.read_html(unicode(electrons_table))[0].iloc[1:]
    df.columns = columns

    return 'electrons', df

def parse_x_gamma_table(x_gamma_table):
    columns = ['type', 'energy', 'intensity', 'dose']

    df = pd.read_html(unicode(x_gamma_table))[0].iloc[1:]
    df.columns = columns

    return 'x_gamma_rays', df

def parse_beta_plus_table(beta_plus_table):
    columns = ['energy', 'end_point_energy', 'intensity', 'dose']

    df = pd.read_html(unicode(beta_plus_table))[0].iloc[1:]
    df.columns = columns

    return 'beta_plus', df


decay_radiation_parsers = {'Electrons': parse_electrons_table,
                           'Beta+': parse_beta_plus_table,
                           'Gamma and X-ray radiation': parse_x_gamma_table}