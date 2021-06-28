import numpy as np


from astropy.io import ascii

import matplotlib as mpl
import matplotlib.cm as cm
mpl.use('Agg')
import matplotlib.pyplot as plt

import flares
import flare.plt as fplt


import flare



cosmo = flare.default_cosmo()

fig, axes = plt.subplots(3, 1, sharex = True, figsize = (3.5, 5))
plt.subplots_adjust(bottom=0.1, left=0.15, top=0.95, wspace=0.0, hspace=0.0)

norm = mpl.colors.Normalize(vmin=8.5, vmax=10.5)

line_styles = ['-.','--','-']


labels = {}
labels['constant10'] = '10\ Myr\ constant\ SF'
labels['constant100'] = '100\ Myr\ constant\ SF'
labels['maxinst'] = 'Max-age\ burst'

for ax, (f1, f2), window in zip(axes, [('Webb.NIRCam.F200W','Webb.NIRCam.F277W'),('Webb.NIRCam.F277W','Webb.NIRCam.F356W'),('Webb.NIRCam.F356W','Webb.NIRCam.F444W')], [(5.6, 6.5),(7.4,8.6),(9.5,10.9)]):

    log10Z = -3
    fesc = 0.0
    for i, sfh_type in enumerate(['maxinst', 'constant10', 'constant100']):

        c = 'k'
        ls = line_styles[i]

        # for fesc, ls in zip([0.0], ['-',':']):

        data = ascii.read(f'data/observed_{sfh_type}_fesc{fesc}_log10Z{log10Z}.dat')
        z = data['z']

        # C = -2.5*np.log10(data[f1]/data[f2])
        C = data[f2]/data[f1]

        ax.plot(z[z<window[0]], C[z<window[0]], c=c, lw=1, ls=ls, alpha=0.2)
        ax.plot(z[(z>window[0])&(z<window[1])], C[(z>window[0])&(z<window[1])], c=c, lw=1, ls=ls, label = rf'$\rm {labels[sfh_type]}$')
        ax.plot(z[z>window[1]], C[z>window[1]], c=c, lw=1, ls=ls, alpha=0.2)
        # ax.plot(data['z'], C, c=c, lw=1, ls=ls, label = rf'$\rm {sfh_type}\ f_{{esc}}={fesc}$')



    log10Z = -3
    fesc = 0.0
    sfh_type = 'constant100'
    for i,log10tau_V in enumerate([-1.0, -0.3, 0.0]):

        c = cm.magma((i+1)/4)
        ls = '-'

        # for fesc, ls in zip([0.0], ['-',':']):

        data = ascii.read(f'data/observed_{sfh_type}_fesc{fesc}_log10Z{log10Z}_log10tauV{log10tau_V}.dat')
        z = data['z']

        # C = -2.5*np.log10(data[f1]/data[f2])
        C = data[f2]/data[f1]

        ax.plot(z[z<window[0]], C[z<window[0]], c=c, lw=1, ls=ls, alpha=0.2)
        ax.plot(z[(z>window[0])&(z<window[1])], C[(z>window[0])&(z<window[1])], c=c, lw=1, ls=ls, label = rf'$\rm \tau_V = {10**log10tau_V:.1f}$')
        ax.plot(z[z>window[1]], C[z>window[1]], c=c, lw=1, ls=ls, alpha=0.2)
        # ax.plot(data['z'], C, c=c, lw=1, ls=ls, label = rf'$\rm {sfh_type}\ f_{{esc}}={fesc}$')



    # no nebular emission
    log10Z = -3
    fesc = 1.0
    sfh_type = 'constant100'

    c = '0.5'
    ls = '-'

    data = ascii.read(f'data/observed_{sfh_type}_fesc{fesc}_log10Z{log10Z}.dat')
    z = data['z']

    # C = -2.5*np.log10(data[f1]/data[f2])
    C = data[f2]/data[f1]

    ax.plot(z[z<window[0]], C[z<window[0]], c=c, lw=1, ls=ls, alpha=0.2)
    ax.plot(z[(z>window[0])&(z<window[1])], C[(z>window[0])&(z<window[1])], c=c, lw=1, ls=ls, label = rf'$\rm f_{{esc}} = 1$')
    ax.plot(z[z>window[1]], C[z>window[1]], c=c, lw=1, ls=ls, alpha=0.2)
    # ax.plot(data['z'], C, c=c, lw=1, ls=ls, label = rf'$\rm {sfh_type}\ f_{{esc}}={fesc}$')



    if ax == axes[0]:ax.legend(fontsize=7, labelspacing = 0.1)

    ax.set_ylim([-0.99, 1.49])
    ax.set_ylim([0, 2.9])
    ax.set_xlim([5., 12])

    if ax == axes[-1]: ax.set_xlabel(r'$\rm z$')
    # ax.set_ylabel(rf'$\rm {f1.split(".")[-1]}-{f2.split(".")[-1]}$')
    ax.set_ylabel(rf'$\rm f_{{ {f1.split(".")[-1]} }}/f_{{ {f2.split(".")[-1]} }}$')

fig.savefig(f'figs/observed_colours.pdf')
fig.clf()
