

# Create a model SED


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import cmasher as cmr

from synthesizer.plt import single
from synthesizer.filters import TopHatFilterCollection
from synthesizer.grid import SpectralGrid
from synthesizer.parametric.sfzh import SFH, ZH, generate_sfzh, generate_instant_sfzh
from synthesizer.parametric.galaxy import SEDGenerator
from unyt import yr, Myr


# -------------------------------------------------
# --- define choise of SPS model and initial mass function (IMF)

model = 'bpass-v2.2.1-bin_chab-300_cloudy-v17.03_log10Uref-2'
grid = SpectralGrid(model)


fig = plt.figure(figsize=(3.5, 3.5))

left = 0.15
height = 0.85
bottom = 0.1
width = 0.8

ax = fig.add_axes((left, bottom, width, height))

a = 0.05
ax.fill_between([3400, 3600], [0, 0], [100, 100], alpha=a, color='r')
ax.fill_between([4150, 4250], [0, 0], [100, 100], alpha=a, color='r')

path_effects = [pe.withStroke(linewidth=2, foreground='white')]

colors = cmr.take_cmap_colors('cmr.rainforest', 4, cmap_range=(0.15, 0.85))


Z = 0.005

Zh = ZH.deltaConstant({'Z': Z})  # constant metallicity


sfh = SFH.Constant({'duration': 100 * Myr})
tauV = 0.0
fesc = 0.0
sfzh = generate_sfzh(grid.log10ages, grid.metallicities, sfh, Zh)

galaxy = SEDGenerator(grid, sfzh)
galaxy.pacman(fesc=fesc, tauV=tauV)  # adds nebular emission

print(galaxy.spectra.keys())

galaxy.spectra['nebular'] = galaxy.spectra['intrinsic']
galaxy.spectra['nebular'].lnu -= galaxy.spectra['stellar'].lnu

for sed_type, ls, alpha, lw in zip(['stellar', 'nebular', 'total'], [':', '-', '-'], [0.4, 0.1, 0.8], [1, 2, 1]):

    sed = galaxy.spectra[sed_type]

    log10bb = sed.get_balmer_break()

    ax.plot(sed.lam, np.log10(sed.lnu)+8, zorder=1, lw=lw,
            c='k', alpha=alpha, ls=ls, label=rf'$\rm {sed_type}$')  # plot SED

    # --- add BB
    l_ = np.array([3500, 4200])
    L_ = np.interp(l_, sed.lam, sed.lnu)*1E8
    ax.plot(l_, np.log10(L_), c='r', alpha=0.3, lw=3)
    ax.scatter(l_, np.log10(L_), c='r', s=5)

    ax.text(3850, np.log10(np.mean(L_)), rf'$\rm {log10bb:.2f}$', color='r',
            fontsize=6, ha='center', rotation=log10bb*80, path_effects=path_effects)


ax.legend(fontsize=8)

ax.set_xlim([2000, 5500])
ax.set_ylim([26.51, 28.99])

ax.legend(fontsize=8)

ax.set_xlabel(r'$\rm \lambda/\regular{\AA} $')
ax.set_ylabel(r'$\rm \log_{10}(L_{\nu}/erg\ s^{-1}\ Hz^{-1}\cdot 10^{8}\ M_{\odot})$')

fig.savefig(f'figs/sed_reprocessing.pdf')

fig.clf()
