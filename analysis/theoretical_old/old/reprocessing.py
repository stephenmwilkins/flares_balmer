import numpy as np

import matplotlib as mpl
import matplotlib.cm as cm
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import cmasher as cmr


from astropy.io import ascii

import flares
import flare.plt as fplt


import flare



fesc = 0.0
log10Z = -3

labels = {}
labels['BB'] = r'$\rm \log_{10}(L_{4200}/L_{3500})$'
labels['beta'] = r'$\rm \beta$'
labels['HbetaEW'] = r'$\rm \log_{10}(H\beta\ EW/\AA)$'


left  = 0.15
height = 0.70
bottom = 0.15
width = 0.8

fig = plt.figure(figsize = (3.5, 3.5))
ax = fig.add_axes((left, bottom, width, height))


colors = cmr.take_cmap_colors('cmr.guppy', 3, cmap_range=(0.3, 0.7))
lw1, lw2 = [1,2]
alpha1, alpha2 = [1, 0.5]
ls1, ls2 = ['-','--']



# --- dust

sfh_type = 'constant'
data = ascii.read(f'data/{sfh_type}_dust.dat')

x = data['log10tau_V']

beta = np.log10(data['FAKE.FAKE.1500']/data['FAKE.FAKE.2500'])/np.log10(1500/2500)-2.0
BB = np.log10(data['FAKE.FAKE.BBb']/data['FAKE.FAKE.BBa'])
HbetaEW = np.log10(data['HbetaEW'])

for k, d, c in zip(['BB', 'beta', 'HbetaEW'], [BB, beta, HbetaEW], colors):
    ax.plot(x, (d - d[0]), c=c, ls=ls1, lw=lw1, alpha = alpha1, label = labels[k])

ax.set_xlim([-2, 1])
ax.set_ylim(-1.5, 1.5)
ax.set_xlabel(r'$\rm log_{10}(\tau_V)$')

# --- fesc

ax2 = ax.twiny()

sfh_type = 'constant'
data = ascii.read(f'data/{sfh_type}_fesc.dat')

x = data['fesc']

beta = np.log10(data['FAKE.FAKE.1500']/data['FAKE.FAKE.2500'])/np.log10(1500/2500)-2.0
BB = np.log10(data['FAKE.FAKE.BBb']/data['FAKE.FAKE.BBa'])
HbetaEW = np.log10(data['HbetaEW'])

for k, d, c in zip(['BB', 'beta', 'HbetaEW'], [BB, beta, HbetaEW], colors):
    ax2.plot(x, (d - d[0]), c=c, ls=ls2, lw=lw2, alpha = alpha2)

ax.axhline(0,lw=3,c='k',alpha=0.1)






ax2.set_xlim([0, 1])
ax2.set_xlabel(r'$\rm f_{esc}$')

ax.set_ylabel(r'$\rm reprocessed - intrinsic$')


handles = [Line2D([0], [0], label = r'$\rm f_{esc}$' , color='k', lw=2, ls='--', alpha=0.5)]
handles += [Line2D([0], [0], label = r'$\rm log_{10}(\tau_V)$' , color='k', lw=1, ls='-', alpha=1.0)]
ax2.add_artist(ax2.legend(loc = 'lower left', handles=handles, fontsize=8, labelspacing = 0.1))


ax.legend(fontsize=8, labelspacing = 0.1)

fig.savefig(f'figs/reprocessing.pdf')
fig.clf()
