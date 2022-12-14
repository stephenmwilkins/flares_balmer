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

# y = 'log10OIIIHbEW'
y_ = 'log10BB'

properties = ['log10sSFR', 'age', 'AFUV', 'log10LV/M*']  # 'log10B','log10M*/LV'
# properties = ['log10sSFR','log10sSFR']

limits[x] = xlimits


limits['log10sSFR'] = [-0.49, 1.99]
limits['log10BB_intrinsic'] = limits['log10BB']
limits['log10LV/M*'] = [18.51, 21.49]

labels['log10BB'] = r'log_{10}(L_{4200}/L_{3500})'
labels['log10BB_intrinsic'] = r'log_{10}(L_{4200}/L_{3500})^{intrinsic}'
labels['log10sSFR'] = r'log_{10}(sSFR_{10}/Gyr^{-1})'
labels['log10LV/M*'] = r'log_{10}[(L_{V}/M_{\star})/erg\ s\ Hz^{-1}\ M_{\odot}^{-1}]'

filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_noparticlesed_v4.hdf5'
flares = analyse.analyse(filename, default_tags=False)


# flares.list_datasets()

# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []

quantities.append({'path': 'Galaxy/Mstar_aperture', 'dataset': f'30',
                  'name': 'Mstar_30', 'log10': True})
quantities.append({'path': f'Galaxy/SFR_aperture/30',
                  'dataset': f'10Myr', 'name': f'SFR_10', 'log10': True})
quantities.append({'path': f'Galaxy/StellarAges',
                  'dataset': 'MassWeightedStellarAge', 'name': 'age', 'log10': True})


quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI',
                  'dataset': 'FUV', 'name': None, 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/Pure_Stellar',
                  'dataset': 'FUV', 'name': 'PSFUV', 'log10': True})

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI',
                  'dataset': 'V', 'name': None, 'log10': True})

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/BB_Wilkins',
                  'dataset': 'DustModelI', 'name': 'BB', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/BB_Wilkins',
                  'dataset': 'Pure_Stellar', 'name': 'BB_intrinsic', 'log10': True})

D = flares.get_datasets(flares.tag_from_zed[z], quantities)
s = D[x] > xlimits[0]

# ----------------------------------------------
# define new quantities
D['log10sSFR'] = D['log10SFR_10']-D['log10Mstar_30']+9
D['AFUV'] = 2.5*(D['log10PSFUV']-D['log10FUV'])
D['log10M*/LV'] = D['log10Mstar_30']-D['log10V']
D['log10LV/M*'] = D['log10V']-D['log10Mstar_30']

for dust_label in ['', '_intrinsic']:

    y = y_+dust_label

    # --- make corner plot
    fig, axes = flares_utility.plt.vlinear(D, y, properties, s, limits=limits, scatter_colour_quantity='log10FUV',
                                           add_correlation_coefficient=True, weighted_median_outline=True, bins=20)

    fig.savefig(f'figs/physical{dust_label}.pdf')
    fig.savefig(f'figs/physical{dust_label}.png')
