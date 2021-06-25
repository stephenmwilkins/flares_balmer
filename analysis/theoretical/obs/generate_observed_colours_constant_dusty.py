


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

fesc = 0.0
sfzh_type = 'constant100'


for log10tau_V in [-1.0, -0.3, 0.0]:

    # --- get SED
    sed = sed_modeller.get_Lnu(sfzh[sfzh_type], {'fesc': fesc, 'log10tau_V': log10tau_V}, dust = ('simple', {'slope':-1}))

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
    ascii.write(data, f'data/observed_{sfzh_type}_fesc{fesc}_log10Z{log10Z}_log10tauV{log10tau_V}.dat', overwrite=True)
