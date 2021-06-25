


import numpy as np

from astropy.io import ascii
from astropy.table import Table, Column, MaskedColumn

import interrogator
from interrogator.sed import models
import interrogator.sed.sfzh

import flare
import flare.filters

cosmo = flare.default_cosmo()


redshifts = np.arange(5., 12., 0.1)

filters = ['Webb.NIRCam.F200W', 'Webb.NIRCam.F277W', 'Webb.NIRCam.F356W', 'Webb.NIRCam.F444W']


# -------------------------------------------------
# --- define choise of SPS model and initial mass function (IMF)
sed_modeller = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')



fesc = 0.0
log10Z = -3


F = flare.filters.add_filters(filters, new_lam = sed_modeller.grid['lam'])


sfzh = {}

# -------------------------------------------------
# --- define star formation and metal enrichment history (sfzh)
log10_duration = 8.
sfzh['constant100'], sfr = interrogator.sed.sfzh.constant(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10_duration': log10_duration, 'log10Z':log10Z, 'log10M*': 1.})

log10_duration = 7.
sfzh['constant10'], sfr = interrogator.sed.sfzh.constant(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10_duration': log10_duration, 'log10Z':log10Z, 'log10M*': 1.})

log10_quenched_duration = 8.
log10_total_duration = 9.
sfzh['quenched'], sfr = interrogator.sed.sfzh.quenched(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10_quenched_duration': log10_quenched_duration, 'log10_total_duration': log10_total_duration,  'log10Z':log10Z, 'log10M*': 1.})



for sfzh_type in ['constant10', 'constant100', 'quenched']:

    for fesc in [0.0,1.0]:

        # --- get SED
        sed = sed_modeller.get_Lnu(sfzh[sfzh_type], {'fesc': fesc}, dust = False)

        d = {f:[] for f in filters}

        for z in redshifts:

            # --- calculate observer frame wavelength
            sed.total.get_fnu(cosmo, z)

            # --- define filters
            F = flare.filters.add_filters(filters, new_lam = sed.total.lamz)

            # --- generates Fnu (broad band fluxes)
            Fnu = sed.total.return_Fnu(F)

            # --- save
            for f in filters:
                d[f].append(Fnu[f])



        # -------------------------------------------------
        # --- save

        data = Table([redshifts]+[d[f] for f in filters], names=['z']+filters)
        ascii.write(data, f'data/observed_{sfzh_type}_fesc{fesc}_log10Z{log10Z}.dat', overwrite=True)
