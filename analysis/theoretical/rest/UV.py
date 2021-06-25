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

fig, ax = fplt.simple()

norm = mpl.colors.Normalize(vmin=8.5, vmax=10.5)


log10Z = -2


for i, sfh_type in enumerate(['constant', 'instant']):

    c = cm.viridis(i/2)

    for fesc, ls in zip([0.0, 1.0], ['-',':']):

        data = ascii.read(f'data/{sfh_type}_fesc{fesc}_log10Z{log10Z}.dat')

        UV = -2.5*np.log10(data['FAKE.FAKE.U']/data['FAKE.FAKE.V'])

        ax.plot(data['log10_duration'], UV, c=c, lw=1, ls=ls, label = f'{sfh_type} f_esc={fesc}')


z = 4
t_uni = np.log10(cosmo.age(z).to('yr').value)
print(t_uni)
ax.axvline(t_uni, label = f'z={z}', c='k', alpha = 0.2)

ax.legend(fontsize=8)

ax.set_ylim([-1., 3])
ax.set_xlim([6., 10])

ax.set_xlabel(r'$\rm log_{10}(t_{SF}/yr)$')
ax.set_ylabel(r'$\rm U-V$')

fig.savefig(f'figs/UV.pdf')
fig.clf()
