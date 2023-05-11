
import numpy as np
import matplotlib.cm as cm
import cmasher as cmr

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

x_limits = [8.01, 11.49]
# x_limits = [28.01, 30.49]

filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_noparticlesed.hdf5'

flares = analyse.analyse(filename, default_tags=False)

# flares.list_datasets()

tags = ['005_z010p000', '006_z009p000', '007_z008p000',
        '008_z007p000', '009_z006p000', '010_z005p000']
zeds = [10, 9, 8, 7, 6, 5]

x = 'log10Mstar_30'
z = 'log10FUV'
y = 'log10BB'


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


flaxes = axes.flatten()

# D4000 = [1.05876971, 0.00664352
# Dsteve = [1.76259398, 0.01156703]


colors = cmr.take_cmap_colors('cmr.infinity', 5, cmap_range=(0.15, 0.85))

x = 10.61
xerr = 0.02
y = np.log10(1.76*(4200**2/3500**2))
print(y)
yerr = np.log10(1.76+0.01156703)-np.log10(1.76)

flaxes[5].errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o',
                   label=r'$\rm Carnall+23$', c=colors[0], markersize=3, elinewidth=1)
z = 7.3
x = 8.7
xerr = 0.1
y = 0.311
yerr = 0.028

flaxes[3].errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o',
                   label=r'$\rm Looser+23$', c=colors[1], markersize=3, elinewidth=1)


z = 10.1
x = 8.1
xerr = 0.3
y = -0.11
yerr = 0.036

flaxes[0].errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o',
                   label=r'$\rm Hsiao+23$', c=colors[2], markersize=3, elinewidth=1)


for ax in axes.flatten():
    ax.legend(fontsize=8)


fig.savefig(f'figs/M.pdf')
fig.savefig(f'figs/M.png')
