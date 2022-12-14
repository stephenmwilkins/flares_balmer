import numpy as np

import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.lines as mlines

# import scipy.stats as stats
from scipy.stats import binned_statistic

import cmasher as cmr

import h5py

from flare.utils import bin_centres
import flare.plt as fplt
import flare.photom as phot
from flare.photom import M_to_lum

from flares_utility.labels import labels
from flares_utility.limits import limits
import flares_utility.plt
import flares_utility.analyse as analyse
import flares_utility.colors

z = 7

x = 'log10FUV'
xlimits = [28.0, 30.49]
diagnostics = ['log10BB', 'log10HbetaEW','beta']
diagnostics = ['log10BB', 'log10OIIIHbEW','beta']
properties = ['log10sSFR', 'log10age', 'log10B', 'AFUV', 'log10M*/LV']
# properties = ['log10sSFR','log10sSFR']

limits[x] = xlimits
# limits['beta'] = [-2.79, -1.61]
# limits['log10HbetaEW'] = [1.26, 2.99]

limits['log10B'] = [-0.99, 1.49]
labels['log10B'] = r'log_{10}(10^{28}\times SFR/L_{FUV})'

filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_noparticlesed.hdf5'
flares = analyse.analyse(filename, default_tags = False)






# flares.list_datasets()

# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []

quantities.append({'path': 'Galaxy/Mstar_aperture', 'dataset': f'30', 'name': 'Mstar_30', 'log10': True})
quantities.append({'path': f'Galaxy/SFR_aperture/30', 'dataset': f'50Myr', 'name': f'SFR_50', 'log10': True})
quantities.append({'path': f'Galaxy/StellarAges', 'dataset': 'MassWeightedStellarAge', 'name': 'age', 'log10': True})
# quantities.append({'path': f'Galaxy/Metallicity', 'dataset': 'MassWeightedGasZ', 'name': None, 'log10': True})

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI', 'dataset': 'FUV', 'name': None, 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/Pure_Stellar', 'dataset': 'FUV', 'name': 'PSFUV', 'log10': True})

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI', 'dataset': 'V', 'name': None, 'log10': True})

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/beta', 'dataset': f'DustModelI', 'name': f'beta', 'log10': False})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/BB_Wilkins', 'dataset': 'DustModelI', 'name': 'BB', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/HI4861', 'dataset': 'EW', 'name': 'HbetaEW', 'log10': True})


quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/HI4861', 'dataset': 'EW', 'name': 'HI4861_EW', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/OIII5007', 'dataset': 'EW', 'name': 'OIII5007_EW', 'log10': False})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/OIII4959', 'dataset': 'EW', 'name': 'OIII4959_EW', 'log10': False})





# --- get all quantities


D = flares.get_datasets(flares.tag_from_zed[z], quantities)
s = D[x]>xlimits[0]

# ----------------------------------------------
# define new quantities
D['log10sSFR'] = D['log10SFR_50']-D['log10Mstar_30']+9
D['AFUV'] = 2.5*(D['log10PSFUV']-D['log10FUV'])
D['log10B'] = D['log10SFR_50']-D['log10FUV']+28
D['log10M*/LV'] = D['log10Mstar_30']-D['log10V']
D['OIIIHbEW'] = D['HI4861_EW']+D['OIII5007_EW']+D['OIII4959_EW']
D['log10OIIIHbEW'] = np.log10(D['OIIIHbEW'])


# --------------------------
# dust attenuated




# --- make corner plot
fig = flares_utility.plt.linear_mcol(D, diagnostics, properties, s, limits = limits, scatter_colour_quantity = 'log10FUV', add_correlation_coefficient = True, weighted_median_outline = True, bins = 20)

fig.savefig(f'figs/physical.pdf')
fig.savefig(f'figs/physical.png')
