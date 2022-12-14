

# Create a model SED


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import cmasher as cmr


import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import cmasher as cmr
import interrogator
from interrogator.sed import models
import interrogator.sed.sfzh
import flare.filters

import flare.plt as fplt
from interrogator.sed.core import rebin



# -------------------------------------------------
# --- define choise of SPS model and initial mass function (IMF)

# SPS = models.SPS('P2/ModSalpeter_100')
SPS = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')
line_modeller = models.lines('BPASSv2.2.1.binary/ModSalpeter_300')



fig = plt.figure(figsize = (3.5,4.5))

left  = 0.15
height = 0.85
bottom = 0.1
width = 0.8

ax = fig.add_axes((left, bottom, width, height))




a = 0.05
ax.fill_between([125, 300], [0,0], [100,100], alpha=a, color='b')
ax.fill_between([340, 360], [0,0], [100,100], alpha=a, color='r')
ax.fill_between([415, 425], [0,0], [100,100], alpha=a, color='r')
for l_ in [500.7,495.9, 486.1]:
    ax.axvline(l_, c='g', alpha = a, lw=2)

# ax.fill_between([220, 280], [0,0], [100,100], alpha=a, color='y')
# ax.fill_between([400, 500], [0,0], [100,100], alpha=a, color='y')


path_effects = [pe.withStroke(linewidth=2, foreground='white')]

colors = cmr.take_cmap_colors('cmr.rainforest', 4, cmap_range=(0.15, 0.85))

for mod, color in zip(['A','B','C','D'], colors):

    # --- young
    if mod == 'A':
        sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': 7., 'log10Z': -3., 'log10M*': 8.})
        dust = False
        params = {'fesc': 0.0}

    # --- average
    if mod == 'B':
        sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': 8., 'log10Z': -3, 'log10M*': 8.})
        dust = False
        params = {'fesc': 0.0}

    # --- average and dusty
    if mod == 'C':
        sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': 8, 'log10Z': -3., 'log10M*': 8.})
        dust = ('simple', {'slope': -1.0})
        params = {'fesc': 0.0, 'log10tau_V': 0.0}

    # --- old
    if mod == 'D':
        sfzh, sfr = interrogator.sed.sfzh.instantaneous(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10age': 9, 'log10Z': -3, 'log10M*': 8.})
        dust = False
        params = {'fesc': 0.0}


    # if mod == 'D':
    #     sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': 7., 'log10Z': -4., 'log10M*': 8.})
    #     dust = False
    #     params = {'fesc': 0.5}

    SED = SPS.get_Lnu(sfzh, params, dust = dust)
    # SED = SPS.get_Lnu(sfzh, {'fesc': 0.0, 'log10tau_V': 0.0}, dust = ('simple', {'slope':-1}))

    # l, L = SED.total.lam, SED.total.lnu
    l, L = rebin(SED.total.lam, SED.total.lnu, 10)

    l /= 10. # convert to nm
    L[l<121.6] = 1E-100 # lyman-alpha break

    ax.plot(l, np.log10(L), zorder = 1, lw=1, c=color, alpha=0.3, label = rf'$\rm {mod}$') # plot SED

    # --- add beta
    # beta = SED.total.return_beta()
    # l_ = np.arange(125, 300, 1)
    # L_ = np.interp(200, l, L)*(l_/200.)**(beta+2.0)
    #
    # ax.plot(l_, np.log10(L_),c='b',alpha=0.3,lw=3)
    # ax.text(180, np.log10(np.interp(200, l, L)), rf'$\rm {beta:.2f}$',color='b',fontsize=6, rotation = (beta+2.0)*22.)


    beta = SED.total.return_beta_spec()
    l_ = np.arange(125, 300, 1)
    L_ = np.interp(200, l, L)*(l_/200.)**(beta+2.0)

    ax.plot(l_, np.log10(L_),c='b',alpha=0.3,lw=3)
    ax.text(180, np.log10(np.interp(200, l, L)), rf'$\rm {beta:.2f}$',color='b',fontsize=6, rotation = (beta+2.0)*22., path_effects=path_effects)


    # --- add BB
    l_ = np.array([350, 420])
    L_ = np.interp(l_, l, L)
    ax.plot(l_, np.log10(L_),c='r',alpha=0.3, lw=3)
    ax.scatter(l_, np.log10(L_),c='r',s=5)

    log10BB = np.log10(L_[1]/L_[0])

    ax.text(385, np.log10(np.mean(L_)), rf'$\rm {log10BB:.2f}$', color='r', fontsize=6, ha = 'center', rotation = log10BB*170, path_effects=path_effects)

    # --- add EW
    EW = 0.
    for line in ['OIII5007','OIII4959','HI4861']:
        q = line_modeller.get_info(sfzh, line, params, dust = dust)
        EW += q['max_EW']*(1-params['fesc'])

    EW /= 10 # convert to nm


    l_ = np.mean([500.7,495.9, 486.1])
    L_ = np.interp(510, l, L)

    ax.plot([l_ - EW/2, l_ + EW/2], [np.log10(L_)]*2, lw=1.5, c='g')
    ax.plot([l_ - EW/2]*2, [np.log10(L_)-0.025, np.log10(L_)+0.025], lw=1.5, c='g')
    ax.plot([l_ + EW/2]*2, [np.log10(L_)-0.025, np.log10(L_)+0.025], lw=1.5, c='g')

    ax.text(495, np.log10(L_)+0.04, rf'$\rm {EW*10:.0f}{{\regular\AA}}$', color='g', fontsize=6, ha = 'center', path_effects=path_effects)


ax.legend(fontsize=8)

ax.set_xlim([100, 600])
ax.set_ylim([25.76, 29.49])

# ax.legend(fontsize=8)

ax.set_xlabel(r'$\rm \lambda/nm $')
ax.set_ylabel(r'$\rm \log_{10}(L_{\nu}/erg\ s^{-1}\ Hz^{-1}\cdot 10^{8}\ M_{\odot})$')

fig.savefig(f'figs/sed.pdf')

fig.clf()
