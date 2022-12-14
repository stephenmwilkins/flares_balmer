

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


fig = plt.figure(figsize=(3.5, 4.5))

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


for mod, color in zip(['A', 'B', 'C', 'D'], colors):

    # --- young
    if mod == 'A':
        sfh = SFH.Constant({'duration': 10 * Myr})
        tauV = 0.0
        fesc = 0.0
        sfzh = generate_sfzh(grid.log10ages, grid.metallicities, sfh, Zh)

    # --- average
    if mod == 'B':
        sfh = SFH.Constant({'duration': 100 * Myr})
        tauV = 0.0
        fesc = 0.0
        sfzh = generate_sfzh(grid.log10ages, grid.metallicities, sfh, Zh)

    # --- average and dusty
    if mod == 'C':
        sfh = SFH.Constant({'duration': 100 * Myr})
        tauV = 1.0
        fesc = 0.0
        sfzh = generate_sfzh(grid.log10ages, grid.metallicities, sfh, Zh)

    # --- old
    if mod == 'D':
        sfzh = generate_instant_sfzh(
            grid.log10ages, grid.metallicities, 9.0, Z)
        tauV = 0.0
        fesc = 0.0

    galaxy = SEDGenerator(grid, sfzh)
    galaxy.pacman(fesc=fesc, tauV=tauV)  # adds nebular emission

    # --- get quanitities

    sed = galaxy.spectra['total']

    log10bb = sed.get_balmer_break()

    # l, L = SED.total.lam, SED.total.lnu
    # l, L = rebin(SED.total.lam, SED.total.lnu, 10)
    #
    # l /= 10.  #  convert to nm
    # L[l < 121.6] = 1E-100  #  lyman-alpha break

    ax.plot(sed.lam, np.log10(sed.lnu)+8, zorder=1, lw=1,
            c=color, alpha=0.3, label=rf'$\rm {mod}$')  # plot SED

    # --- add BB
    l_ = np.array([3500, 4200])
    L_ = np.interp(l_, sed.lam, sed.lnu)*1E8
    ax.plot(l_, np.log10(L_), c='r', alpha=0.3, lw=3)
    ax.scatter(l_, np.log10(L_), c='r', s=5)

    ax.text(3850, np.log10(np.mean(L_)), rf'$\rm {log10bb:.2f}$', color='r',
            fontsize=6, ha='center', rotation=log10bb*80, path_effects=path_effects)


ax.legend(fontsize=8)

ax.set_xlim([2000, 5500])
ax.set_ylim([25.76, 29.49])

# ax.legend(fontsize=8)

ax.set_xlabel(r'$\rm \lambda/\regular{\AA} $')
ax.set_ylabel(r'$\rm \log_{10}(L_{\nu}/erg\ s^{-1}\ Hz^{-1}\cdot 10^{8}\ M_{\odot})$')

fig.savefig(f'figs/sed.pdf')

fig.clf()
