
import matplotlib as mpl
import matplotlib.cm as cm
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import numpy as np

import pickle

import flares
import flares_analysis
import flare.plt as fplt



O = pickle.load(open('observed_colours_8.5.p','rb'))
Ohm = pickle.load(open('observed_colours_9.5.p','rb'))
Ohm2 = pickle.load(open('observed_colours_10.5.p','rb'))


filters = [f'Webb.NIRCAM.{f}' for f in ['F200W','F277W','F356W','F444W']]


# --- define colour scale
norm = mpl.colors.Normalize(vmin=5, vmax=10)


fig, axes = plt.subplots(len(filters)-1, 1, figsize = (3.5,3), sharex = True)
plt.subplots_adjust(left=0.15, top=0.95, bottom=0.15, right=0.9, wspace=0.0, hspace=0.0)


for f1, f2, ax in zip(filters[:-1], filters[1:], axes):

    for tag, z in zip(O.keys(), np.arange(10, 4, -1)):

        c = cm.plasma(norm(z))

        ax.fill_between(O[tag]['z'], O[tag][f'{f1}_{f2}_P16'], O[tag][f'{f1}_{f2}_P84'], color=c, alpha=0.2)
        ax.plot(O[tag]['z'], O[tag][f'{f1}_{f2}_P50'], c=c, lw=1)

        ax.plot(Ohm[tag]['z'], Ohm[tag][f'{f1}_{f2}_P50'], c=c, lw=1, ls='-.')
        if len(Ohm2[tag][f'{f1}_{f2}_P50']) == len(Ohm2[tag]['z']):
            ax.plot(Ohm2[tag]['z'], Ohm2[tag][f'{f1}_{f2}_P50'], c=c, lw=1, ls=':')

    ax.set_xlim(4.5, 10.5)
    ax.set_ylim(0.6, 2.6)
    ax.set_ylabel(rf"$\rm f_{{ {f2.split('.')[-1]} }}/f_{{ {f1.split('.')[-1]} }} $", fontsize=8)

axes[-1].set_xlabel(rf'$\rm z $')


observations = []
observations.append({'id': 'MACS1149-JD1', 'ref': 'Hashimoto', 'z': 9.1096, 'F444W/F356W': [2.3, [1.9, 2.7]]})

for obs in observations:

    axes[-1].scatter(obs['z'], obs['F444W/F356W'][0], c='k')
    axes[-1].plot([obs['z']]*2, obs['F444W/F356W'][1], c='k', lw=1)




handles = [Line2D([0], [0], label = r'$\rm log_{10}(M_{\star}/M_{\odot})>8.5$' , color='k', lw=1, ls='-')]
handles += [Line2D([0], [0], label = r'$\rm log_{10}(M_{\star}/M_{\odot})>9.5$' , color='k', lw=1, ls='-.')]
handles += [Line2D([0], [0], label = r'$\rm log_{10}(M_{\star}/M_{\odot})>10.5$' , color='k', lw=1, ls=':')]


axes[0].legend(handles=handles, fontsize=7, labelspacing = 0.1)

fig.savefig(f'figs/colour_evolution.pdf')
