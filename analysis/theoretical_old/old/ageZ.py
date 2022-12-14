import numpy as np


from astropy.io import ascii

import matplotlib as mpl
import matplotlib.cm as cm
import cmasher as cmr
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import flares
import flare.plt as fplt


import flare



fesc = 0.0
log10Z = -3
tau_V = 0.0

fig, axes = plt.subplots(3, 1, figsize = (3.5,5), sharex = True)
plt.subplots_adjust(left=0.15, top=0.975, bottom=0.1, right=0.95, wspace=0.02, hspace=0.02)


c = '0.5'

log10Zs = [-3, -2]
colors = cmr.take_cmap_colors('cmr.watermelon', len(log10Zs), cmap_range=(0.25, 0.75))

line_styles = ['-','--','-.',':']


for log10Z, c in zip(log10Zs, colors):

    for i, sfh_type in enumerate(['exp100', 'constant', 'exp-100', 'instant']):

        ls = line_styles[i]

        data = ascii.read(f'data/{sfh_type}_fesc{fesc}_log10Z{log10Z}_tauV{tau_V}.dat')
        if sfh_type == 'instant':
            x = data['log10_duration']
        elif sfh_type == 'constant':
            x = data['log10_duration']-np.log10(2)
        elif sfh_type == 'exp100':
            tau = 100
            t = 10**(data['log10_duration']-6.)/tau
            x = -np.log(0.5*(1+np.exp(-t)))
            x = np.log10(x) + np.log10(tau) + 6
        elif sfh_type == 'exp-100':
            tau = -100
            t = 10**(data['log10_duration']-6.)/np.fabs(tau)
            x = np.log(0.5*(1+np.exp(t)))
            x = np.log10(x) + np.log10(np.fabs(tau)) + 6


        beta = np.log10(data['FAKE.FAKE.1500']/data['FAKE.FAKE.2500'])/np.log10(1500/2500)-2.0
        BB = np.log10(data['FAKE.FAKE.BBb']/data['FAKE.FAKE.BBa'])
        HbetaEW = np.log10(data['HbetaEW'])

        for ax, d in zip(axes, [BB, beta, HbetaEW]):
            ax.plot(x, d, c=c, lw=1, ls=ls)




c='0.5'
handles = [Line2D([0], [0], label = fr'$\rm instantaneous\ burst$' , color=c, lw=1, ls=':')]
handles += [Line2D([0], [0], label = fr'$\rm exp\ decreasing\ SF\ (\tau=-100\ Myr)$' , color=c, lw=1, ls='-.')]
handles += [Line2D([0], [0], label = fr'$\rm constant\ SF$' , color=c, lw=1, ls='--')]
handles += [Line2D([0], [0], label = fr'$\rm exp\ increasing\ SF\ (\tau=100\ Myr)$' , color=c, lw=1, ls='-')]

handles += [Line2D([0], [0], label = fr'$\rm \log_{{10}}(Z)={log10Z}$' , color=c, lw=1, ls='-') for log10Z, c in zip(log10Zs, colors)]

leg = axes[0].legend(loc = 'upper left', handles=handles, fontsize=8, labelspacing = 0.1)
axes[0].add_artist(leg)


limits = {}
limits['BB'] = [-0.3, 0.9]
limits['beta'] = [-3.5, 0.9]
limits['HbetaEW'] = [2.9, 0.1]

labels = {}
labels['BB'] = r'$\rm \log_{10}(L_{4200}/L_{3500})$'
labels['beta'] = r'$\rm \beta$'
labels['HbetaEW'] = r'$\rm \log_{10}(H\beta\ EW/\AA)$'

for ax, d in zip(axes, ['BB','beta','HbetaEW']):
    ax.set_ylim(limits[d])
    ax.set_xlim([6., 9.0])
    ax.set_ylabel(labels[d])

ax.set_xlabel(r'$\rm log_{10}(age/yr)$')


fig.savefig(f'figs/ageZ.pdf')
fig.clf()
