
import numpy as np
import matplotlib.cm as cm

import scipy.stats as stats


import cmasher as cmr

import h5py

import flare.plt as fplt
import flare.photom as phot
from flare.photom import M_to_lum
import flares_utility.limits
import flares_utility.plt
import flares_utility.analyse as analyse
import flares_utility.stats

# x_limits = [8.5, 11.5]
x_limits = [28.01, 30.49]

filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_noparticlesed.hdf5'

flares = analyse.analyse(filename, default_tags=False)

# flares.list_datasets()

tags = ['005_z010p000', '006_z009p000', '007_z008p000',
        '008_z007p000', '009_z006p000', '010_z005p000']
zeds = [10, 9, 8, 7, 6, 5]


x = 'log10Mstar_30'
x = 'log10FUV'
y = 'log10BB'
z = 'log10Mstar_30'

# flares.list_datasets()

# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []

quantities.append({'path': 'Galaxy/Mstar_aperture', 'dataset': f'30',
                  'name': 'Mstar_30', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI',
                  'dataset': 'FUV', 'name': None, 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/BB_Wilkins',
                  'dataset': 'DustModelI', 'name': 'BB', 'log10': True})


D = {}
s = {}
s = {}

for tag, redshift in zip(tags, zeds):

    # --- get quantities (and weights and deltas)
    D[redshift] = flares.get_datasets(tag, quantities)
    s[redshift] = D[redshift][x] > x_limits[0]


limits = flares_utility.limits.limits

# limits[y] = [10, 4000]
limits[y] = [-0.3, 0.6]
limits[x] = x_limits


# cmap = cmr.get_sub_cmap('cmr.sapphire', 0.4, 1.0)
cmap = cmr.bubblegum

# fig, axes = flares_utility.plt.linear_redshift(
#     D, zeds, x, y, s, limits=limits, scatter_colour_quantity=z, scatter_cmap=cmap, bins=15, rows=2, lowz=True, add_weighted_range=False)

fig, axes = flares_utility.plt.linear_redshift(
    D, zeds, x, y, s, limits=limits, scatter=False, bins=15, rows=2, lowz=True, add_weighted_range=True)

# for ax in axes.flatten():
#     ax.set_yscale('log')


fig.savefig(f'figs/L.pdf')
fig.savefig(f'figs/L.png')
