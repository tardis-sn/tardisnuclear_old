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
        print "Reading",
        for nuclear_name in isotope_list:
            print nuclear_name
            isotope_nuclear_data = get_decay_radiation(nuclear_name)
            for data_name, data_table in isotope_nuclear_data.items():
                if (('energy' in data_table.columns)
                    and ('intensity' in data_table.columns)):
                    data_table['energy_per_decay'] = (data_table.energy *
                                                      data_table.intensity)
                    isotope_nuclear_data[
                        ('total_{0}_energy_per_decay'.format(data_name))] = \
                        data_table.energy_per_decay.sum()

                    leptons_energy = (
                        isotope_nuclear_data.get(
                            'total_beta_plus_energy_per_decay', 0.0) +
                        isotope_nuclear_data.get(
                            'total_beta_minus_energy_per_decay', 0.0) +
                        isotope_nuclear_data.get(
                            'total_electrons_energy_per_decay', 0.0))
                    isotope_nuclear_data['total_lepton_energy_per_decay'] = (
                        leptons_energy
                    )


            decay_radiation[nuclear_name] = isotope_nuclear_data

        return decay_radiation

