


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
# redshifts = np.arange(5., 12., 1)

filters = ['Webb.NIRCam.F200W', 'Webb.NIRCam.F277W', 'Webb.NIRCam.F356W', 'Webb.NIRCam.F444W']


# -------------------------------------------------
# --- define choise of SPS model and initial mass function (IMF)
sed_modeller = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')


log10Z = -3
fesc = 0.0

ages = []
d = {f:[] for f in filters}

for z in redshifts:

    # --- determine age of Universe

    age = cosmo.age(z).to('yr').value
    ages.append(age)
    print(z, age)

    sfzh, sfr = interrogator.sed.sfzh.instantaneous(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10age': np.log10(age), 'log10Z':log10Z, 'log10M*': 1.})

    # --- get SED
    sed = sed_modeller.get_Lnu(sfzh, {'fesc': fesc}, dust = False)

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
ascii.write(data, f'data/observed_maxinst_fesc{fesc}_log10Z{log10Z}.dat', overwrite=True)
