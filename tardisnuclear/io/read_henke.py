import requests
from bs4 import BeautifulSoup

base_url = 'http://henke.lbl.gov/cgi-bin/pert_cgi.pl'


def get_photo_absorption_cross_section(element_code, energy):
    data = requests.get(base_url, data=dict(
        Element=element_code, Energy=energy))
    bs = BeautifulSoup(data.text)
    for item in bs.find_all('li'):
        if item.text.startswith('Photo'):
            return float(item.text.split(':')[1].strip().replace('cm^2/g', ''))


