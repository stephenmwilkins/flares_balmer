

# Create a model SED


import numpy as np
import matplotlib.pyplot as plt

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


# -------------------------------------------------
# --- define star formation and metal enrichment history (sfzh)

sfzh, sfr = interrogator.sed.sfzh.constant(SPS.grid['log10age'], SPS.grid['log10Z'] , {'log10_duration': 9., 'log10Z': -3., 'log10M*': 8.})

# plt.imshow(sfzh)
# plt.show()


print('star formation rate: {0}'.format(sfr))

SED = SPS.get_Lnu(sfzh, {'fesc': 0.0}, dust = False)
# SED = SPS.get_Lnu(sfzh, {'fesc': 0.0, 'log10tau_V': 0.0}, dust = ('simple', {'slope':-1}))


fig, ax = fplt.simple()

# l, L = SED.total.lam, SED.total.lnu
l, L = rebin(SED.total.lam, SED.total.lnu, 10)

l /= 10. # convert to nm


L[l<121.6] = 1E-100 # lyman-alpha break

log10L = np.log10(L)

print(np.min(log10L), np.max(log10L))



# add beta
beta = SED.total.return_beta()
l_ = np.array([125, 300])
L_ = np.interp(200, l, L)*(l_/200.)**(beta+2.0)

ax.plot(np.log10(l_), np.log10(L_),c='b',alpha=0.2,lw=4)



# add BB
l_ = np.array([350, 450])
L_ = np.interp(l_, l, L)
ax.plot(np.log10(l_), np.log10(L_),c='r',alpha=0.3, lw=2)

# --- This work
ax.fill_between([340, 360], [0,0], [100,100], alpha=0.1, color='r')
ax.fill_between([415, 425], [0,0], [100,100], alpha=0.1, color='r')

for l_ in [500.7,495.9, 486.1]:
    ax.axvline(np.log10(l_), c='g', alpha = 0.2, lw=2)






ax.plot(np.log10(l), np.log10(L), zorder = 1, lw=1, c='0.5') # plot SED



ax.set_xlim([2., 2.79])
ax.set_ylim([26.51, 28.99])

# ax.legend(fontsize=8)

ax.set_xlabel(r'$\rm \log_{10}(\lambda/nm) $')
ax.set_ylabel(r'$\rm \log_{10}(L_{\nu}/erg\ s^{-1}\ Hz^{-1}\cdot 10^{8}\ M_{\odot})$')

fig.savefig(f'figs/sed.pdf')

fig.clf()
