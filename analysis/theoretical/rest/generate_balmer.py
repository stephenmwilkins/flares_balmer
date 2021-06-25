


import numpy as np

from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn

import interrogator
from interrogator.sed import models
import interrogator.sed.sfzh

import flare.filters




# -------------------------------------------------
# --- define choise of SPS model and initial mass function (IMF)
sed_modeller = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')
line_modeller = models.lines('BPASSv2.2.1.binary/ModSalpeter_300')



log10_durations = np.arange(6., 10.1, 0.1)
filters = ['FAKE.FAKE.1500', 'FAKE.FAKE.2500', 'FAKE.FAKE.BBa', 'FAKE.FAKE.BBb', 'FAKE.FAKE.D4000a', 'FAKE.FAKE.D4000b', 'FAKE.FAKE.U', 'FAKE.FAKE.V']
F = flare.filters.add_filters(filters, new_lam = sed_modeller.grid['lam'])



# for sfh_type in ['exp100', 'constant', 'instant']:
for sfh_type in ['instant']:

    for log10Z in [-3, -2]:

        for fesc in [0.0, 0.5, 0.9, 1.0]:

            for log10tau_V in [-5., -1.0, -0.3, 0.0]:

                tau_V = np.round(10**log10tau_V, 1)

                print(sfh_type, log10Z, fesc, tau_V)

                d = {}
                for f in filters: d[f] = []
                d['HbetaEW'] = []

                for log10_duration in log10_durations:

                    # -------------------------------------------------
                    # --- define star formation and metal enrichment history (sfzh)

                    if sfh_type == 'constant':
                        sfzh, sfr = interrogator.sed.sfzh.constant(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10_duration': log10_duration, 'log10Z':log10Z, 'log10M*': 1.})

                    if sfh_type == 'instant':
                        sfzh, sfr = interrogator.sed.sfzh.instantaneous(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10age': log10_duration, 'log10Z':log10Z, 'log10M*': 1.})

                    if sfh_type == 'exp100':
                        sfzh, sfr = interrogator.sed.sfzh.exponential(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10_duration': log10_duration, 'tau_Myr': 100., 'log10Z':log10Z, 'log10M*': 1.})


                    # --- get SED
                    # sed = sed_modeller.get_Lnu(sfzh, {'fesc': fesc}, dust = False)
                    sed = sed_modeller.get_Lnu(sfzh, {'fesc': fesc, 'log10tau_V': log10tau_V}, dust = ('simple', {'slope':-1}))

                    # --- get FUV and NUV luminosities
                    L = sed.total.return_Lnu(F)

                    for f in filters:
                        d[f].append(L[f])

                    # --- get emission line quantities, specifically Hbeta
                    Hbeta = line_modeller.get_info(sfzh, 'HI4861', {'fesc': fesc, 'log10tau_V': log10tau_V}, dust = ('simple', {'slope':-1}))

                    d['HbetaEW'].append(Hbeta['EW'])




                # -------------------------------------------------
                # --- save

                data = Table([log10_durations]+[d[f] for f in filters]+[d['HbetaEW']], names=['log10_duration']+filters+['HbetaEW'])
                ascii.write(data, f'data/{sfh_type}_fesc{fesc}_log10Z{log10Z}_tauV{tau_V}.dat', overwrite=True)
