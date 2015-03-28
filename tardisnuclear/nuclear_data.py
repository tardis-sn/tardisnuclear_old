from tardisnuclear.io import get_decay_radiation
from pyne import nucname

class NuclearData(object):

    def __init__(self, isotope_list):
        self.data = self._get_decay_radiation_data(isotope_list)


    def __getitem__(self, item):
        try:
            isotope = item.title()
        except AttributeError:
            try:
                nucname.name(item)
            except RuntimeError:
                raise ValueError('item is neither a integer or string that '
                                 'identifies an isotope')
            else:
                isotope = nucname.name(item)

        return self.data[isotope]

    @staticmethod
    def _get_decay_radiation_data(isotope_list):
        decay_radiation = {}
        for nuclear_name in isotope_list:
            decay_radiation[nuclear_name] = get_decay_radiation(nuclear_name)
        return decay_radiation

