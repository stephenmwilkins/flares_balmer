


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



fesc = 0.0
log10Z = -3
log10_duration = 8.
log10tau_V = -5.

fescs = np.arange(0., 1, 0.01)

for sfh_type in ['exp-100', 'exp100', 'constant', 'instant']:

    print(sfh_type)

    # --- output files
    d = {}
    for f in filters: d[f] = []
    d['HbetaEW'] = []


    # -------------------------------------------------
    # --- define star formation and metal enrichment history (sfzh)

    if sfh_type == 'constant':
        sfzh, sfr = interrogator.sed.sfzh.constant(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10_duration': log10_duration, 'log10Z':log10Z, 'log10M*': 1.})

    if sfh_type == 'instant':
        sfzh, sfr = interrogator.sed.sfzh.instantaneous(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10age': log10_duration, 'log10Z':log10Z, 'log10M*': 1.})

    if sfh_type == 'exp100':
        sfzh, sfr = interrogator.sed.sfzh.exponential(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10_duration': log10_duration, 'tau_Myr': 100., 'log10Z':log10Z, 'log10M*': 1.})

    if sfh_type == 'exp-100':
        sfzh, sfr = interrogator.sed.sfzh.exponential(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10_duration': log10_duration, 'tau_Myr': -100., 'log10Z':log10Z, 'log10M*': 1.})

    d['Hbeta_EW'] = []
    d['OIIIr_EW'] = []
    d['OIIIb_EW'] = []
    d['OIIIHbeta_EW'] = []


    for fesc in fescs:

        print(fesc)

        # --- get SED
        # sed = sed_modeller.get_Lnu(sfzh, {'fesc': fesc}, dust = False)
        sed = sed_modeller.get_Lnu(sfzh, {'fesc': fesc, 'log10tau_V': log10tau_V}, dust = ('simple', {'slope':-1}))

        # --- get FUV and NUV luminosities
        L = sed.total.return_Lnu(F)

        for f in filters:
            d[f].append(L[f])

        OIIIr = line_modeller.get_info(sfzh, 'OIII5007', {'fesc': fesc, 'log10tau_V': log10tau_V}, dust = ('simple', {'slope':-1}))
        OIIIb = line_modeller.get_info(sfzh, 'OIII4959', {'fesc': fesc, 'log10tau_V': log10tau_V}, dust = ('simple', {'slope':-1}))
        Hbeta = line_modeller.get_info(sfzh, 'HI4861', {'fesc': fesc, 'log10tau_V': log10tau_V}, dust = ('simple', {'slope':-1}))

        d['Hbeta_EW'].append(Hbeta['EW'])
        d['OIIIr_EW'].append(OIIIr['EW'])
        d['OIIIb_EW'].append(OIIIb['EW'])
        d['OIIIHbeta_EW'].append(OIIIb['EW']+OIIIr['EW']+Hbeta['EW'])

    # -------------------------------------------------
    # --- save


    data = Table([fescs]+[d[f] for f in filters]+[d[l] for l in ['Hbeta_EW','OIIIr_EW','OIIIb_EW','OIIIHbeta_EW']], names=['fesc']+filters+['Hbeta_EW','OIIIr_EW','OIIIb_EW','OIIIHbeta_EW'])
    ascii.write(data, f'data/{sfh_type}_fesc.dat', overwrite=True)
