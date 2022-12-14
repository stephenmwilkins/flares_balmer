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

x = 'log10FUV'
property = 'log10BB'
xlimits = [28.0, 30.49]
ylimits = [-0.05, 0.35]


filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_noparticlesed.hdf5'
flares = analyse.analyse(filename, default_tags=False)

# redshifts = flares.zeds
redshifts = [5, 7, 9]
colors = cmr.take_cmap_colors('cmr.infinity', 3, cmap_range=(0.15, 0.85))


# flares.list_datasets()

# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI',
                  'dataset': 'FUV', 'name': None, 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/Pure_Stellar',
                  'dataset': 'FUV', 'name': 'PSFUV', 'log10': True})

for t in ['DustModelI', 'Intrinsic', 'Pure_Stellar']:
    quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/BB_Wilkins',
                      'dataset': t, 'name': f'BB_{t}', 'log10': True})


# --- get all quantities
D = {}
s = {}
for z in redshifts:
    D[z] = flares.get_datasets(flares.tag_from_zed[z], quantities)
    s[z] = D[z][x] > xlimits[0]


bin_edges = np.arange(*xlimits, 0.25)


fig, ax = fplt.simple()


ptype_ = ['Pure_Stellar', 'Intrinsic', 'DustModelI']
ls_ = [':', '--', '-']
labels_ = [r'$\rm stellar$', r'$\rm stellar\ +\ nebular$',
           r'$\rm stellar\ +\ nebular\ +\ dust$']

for z, c in zip([5, 7, 9], colors):
    for t, ls in zip(ptype_, ls_):

        Y, bin_edges, _ = binned_statistic(
            D[z]['log10FUV'][s[z]], D[z][f'{property}_{t}'][s[z]], statistic='median', bins=bin_edges)

        ax.plot(bin_centres(bin_edges), Y, c=c, lw=1, ls=ls)

handles = [mlines.Line2D([], [], color='0.5', ls=ls, lw=1, label=label)
           for ls, label in zip(ls_, labels_)]
handles += [mlines.Line2D([], [], color=c, ls='-', lw=1,
                          label=rf'$\rm z={z:.0f}$') for z, c in zip(redshifts, colors)]

ax.legend(handles=handles, fontsize=7, labelspacing=0.0)

# ax.xaxis.set_major_locator(plt.MaxNLocator(6))
ax.set_xticks(np.arange(28, 31, 1.0))
ax.set_ylim(ylimits)
ax.set_yticks([0.0, 0.1, 0.2, 0.3])

ax.set_xlabel(r'$\rm log_{10}(L_{FUV}/erg\ s^{-1}\ Hz^{-1})$')
ax.set_ylabel(rf'$\rm {labels[property]}$')
fig.savefig(f'figs/reprocessing.pdf')
fig.clf()
