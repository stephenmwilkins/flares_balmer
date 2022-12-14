

# Import CMasher to register colormaps
import cmasher as cmr

import numpy as np
import matplotlib.cm as cm

import pickle

from astropy.io import ascii
import flares_utility.analyse as analyse
import flares_utility.limits
import flares_utility.labels
import flares_utility.plt

from flare.obs.literature.general import get_binned
import flare.obs.literature.OIIIHb as OIIIHb_observations


filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_noparticlesed.hdf5'
flares = analyse.analyse(filename, default_tags = False)

x_limit = 28.0
xlimits = [x_limit, 30.49]


x = 'log10FUV'
y = 'OIIIHbEW'

cmap = cmr.infinity

limits = flares_utility.limits.limits
limits[x][0] = x_limit
limits['log10sSFR'] = [-0.99, 1.49]


labels = flares_utility.labels.labels


# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []

quantities.append({'path': 'Galaxy/Mstar_aperture', 'dataset': f'30', 'name': 'Mstar_30', 'log10': True})
quantities.append({'path': f'Galaxy/SFR_aperture/30', 'dataset': f'50Myr', 'name': f'SFR_50', 'log10': True})

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI', 'dataset': 'FUV', 'name': None, 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/HI4861', 'dataset': 'EW', 'name': 'HI4861_EW', 'log10': False})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/OIII5007', 'dataset': 'EW', 'name': 'OIII5007_EW', 'log10': False})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/OIII4959', 'dataset': 'EW', 'name': 'OIII4959_EW', 'log10': False})

redshifts = [7,8]
observations = [OIIIHb_observations.observed[redshift] for redshift in redshifts]
colors = cmr.take_cmap_colors('cmr.voltage', len(observations), cmap_range=(0.15, 0.85))
i = 0

for redshift in redshifts:

    # --- get quantities (and weights and deltas)
    D = flares.get_datasets(flares.tag_from_zed[redshift], quantities)

    D['OIIIHbEW'] = D['HI4861_EW']+D['OIII5007_EW']+D['OIII4959_EW']
    D['log10OIIIHbEW'] = np.log10(D['OIIIHbEW'])
    D['log10sSFR'] = D['log10SFR_50'] - D['log10Mstar_30'] + 9

    s = D['log10FUV']>x_limit

    fig, ax, hax = flares_utility.plt.simple_whist(D, x, y, s, limits = limits, labels = labels)

    # --- add observational comparisons

    markers = ['o','d','s']


    x_ = 'log10L1500'
    y_ = 'OIIIHbEW'

    if redshift in OIIIHb_observations.observed.keys():

        for obs, marker in zip(OIIIHb_observations.observed[redshift], markers):

            c = colors[i]
            i += 1

            if obs.dt == 'io':
                ax.scatter(obs.i[x_], obs.i[y_], c=c, s=3, marker=marker, label = rf'$\rm {obs.label}\ (individual)$')
                hax.hist(obs.i[y_], bins=30, orientation='horizontal',color=c,density=True, alpha = 0.5,lw=0)

                b = get_binned(obs.i, x_, y_, redshift, np.arange(27, 31, 0.5))

                ax.errorbar(b[x_], b[y_], xerr = b[x_+'_err'], yerr = b[y_+'_err'], fmt='o', elinewidth=1, ms=5, label = rf'$\rm {obs.label}\ (binned)$', c='k', mec='k', mfc=c, zorder = 4)


    ax.legend(fontsize=7, handletextpad = 0.0)
    ax.set_xticks(np.arange(29, 31, 1.0))
    ax.text(1.1, 0.9,rf'$\rm z={redshift}$', fontsize = 9, transform=ax.transAxes)

    fig.savefig(f'figs/OIIIHb_observations_{redshift}.pdf')
