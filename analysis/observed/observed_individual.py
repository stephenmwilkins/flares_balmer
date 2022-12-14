
import numpy as np
import matplotlib.cm as cm
import matplotlib.patheffects as pe
import cmasher as cmr

import flare.plt as fplt
from flare.utils import bin_centres
import flares_utility.analyse as analyse
import flares_utility.limits
import flares_utility.plt

from flare.obs.literature.general import get_binned
import flare.obs.literature.beta as beta_observations
import flare.obs.literature.OIIIHb as OIIIHb_observations

# ----------------------------------------------------------------------
# --- open data


filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_noparticlesed.hdf5'

flares = analyse.analyse(filename, default_tags = False)

x = 'log10FUV'
xlimits = [28.0, 30.49]
properties = ['log10BB', 'log10HI4861_EW','beta','log10OIIIHbEW']
# properties = ['beta']

z = 'log10Mstar_30'


# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []

quantities.append({'path': 'Galaxy/Mstar_aperture', 'dataset': f'30', 'name': 'Mstar_30', 'log10': True})

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/beta', 'dataset': f'DustModelI', 'name': f'beta', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/BB_Wilkins', 'dataset': 'DustModelI', 'name': 'BB', 'log10': True})

# quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/HI4861', 'dataset': 'EW', 'name': 'HbetaEW', 'log10': True})

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/HI4861', 'dataset': 'EW', 'name': 'HI4861_EW', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/OIII5007', 'dataset': 'EW', 'name': 'OIII5007_EW', 'log10': False})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/OIII4959', 'dataset': 'EW', 'name': 'OIII4959_EW', 'log10': False})


quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI', 'dataset': 'FUV', 'name': None, 'log10': True})







D = {}
s = {}

for tag, redshift in zip(flares.tags, flares.zeds):

    # --- get quantities (and weights and deltas)
    D[redshift] = flares.get_datasets(tag, quantities)

    D[redshift]['OIIIHbEW'] = D[redshift]['HI4861_EW']+D[redshift]['OIII5007_EW']+D[redshift]['OIII4959_EW']
    D[redshift]['log10OIIIHbEW'] = np.log10(D[redshift]['OIIIHbEW'])

    s[redshift] = D[redshift][x]>xlimits[0]


limits = flares_utility.limits.limits

limits[x] = xlimits


for y in properties:

    fig, axes = flares_utility.plt.linear_redshift_dual_wmag(D, flares.zeds, x, y, z, s, limits = limits)

    if y == 'beta':
        observations = beta_observations
        x_ = 'log10LFUV'
        y_ = 'beta'
    elif y == 'log10OIIIHbEW':
        observations = OIIIHb_observations
        x_ = 'log10L1500'
        y_ = 'log10OIIIHbEW'
    else:
        observations = None

    if observations:

        observation_list = []
        added_to_legend = []

        for redshift in flares.zeds:
            if redshift in observations.observed.keys():
                for obs in observations.observed[redshift]:
                    observation_list.append(obs.label)

        observation_list = set(observation_list)

        cmap = 'cmr.voltage'
        colors = dict(zip(observation_list, cmr.take_cmap_colors(cmap, len(observation_list), cmap_range=(0.15, 0.85)))) #, cmap_range=(0.15, 0.85)
        markers = dict(zip(observation_list,  ['o','v','D','s','^','p','h','d']))

        for i, redshift in enumerate(flares.zeds):
            ax = axes[1,i]
            if redshift in observations.observed.keys():
                for obs in observations.observed[redshift]:

                    id = obs.label

                    if id not in added_to_legend:
                        label = rf'$\rm {id}$'
                        added_to_legend.append(id)
                    else:
                        label = None

                    print(redshift, label)

                    # if obs.dt == 'io':
                    #     ax.scatter(obs.d[x_], obs.d[y_], c=[colors[id]], s=3, marker=markers[id], label = label)

                    if obs.dt == 'io':
                        d = get_binned(obs.i, x_, y_, redshift, np.arange(27, 31, 0.5))
                    if obs.dt == 'binned':
                        d = obs.d[redshift]

                    if type(d[x_])==np.ndarray:
                        if y_+'_err' in d.keys():
                            ax.errorbar(d[x_], d[y_], xerr = d[x_+'_err'], yerr = d[y_+'_err'], fmt='o', elinewidth=1, ms=3, label = label, c=colors[id], mec=colors[id], mfc=colors[id], zorder = 4)
                        else:
                            ax.scatter(d[x_], d[y_], label = label, zorder = 2)







                ax.legend(fontsize=7, handletextpad = 0.0, loc = 'upper right')

    fig.savefig(f'figs/observed_{y}.pdf')
    fig.savefig(f'figs/observed_{y}.png')
