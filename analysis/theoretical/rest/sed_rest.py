




import numpy as np
import matplotlib.pyplot as plt

import flare.plt as fplt

import interrogator
from interrogator.sed import models
import interrogator.sed.sfzh
import flare.filters




# -------------------------------------------------
# --- define choise of SPS model and initial mass function (IMF)

sed_modeller = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')







fig, ax = fplt.simple()


# --- constant, no dust
sfzh, sfr = interrogator.sed.sfzh.constant(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10_duration': 8., 'log10Z': -3., 'log10M*': 0.})
sed = sed_modeller.get_Lnu(sfzh, {'fesc': 0.0}, dust = False)
ax.plot(sed.total.lam, np.log10(sed.total.lnu), c='0.5', lw=1)

# --- constant A_V = 1
sfzh, sfr = interrogator.sed.sfzh.constant(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10_duration': 8., 'log10Z': -3., 'log10M*': 0.})
sed = sed_modeller.get_Lnu(sfzh, {'fesc': 0.0, 'log10tau_V': 0.0}, dust = ('simple', {'slope':-1}))
ax.plot(sed.total.lam, np.log10(sed.total.lnu), c='0.5', lw=1)

# --- quenched
sfzh, sfr = interrogator.sed.sfzh.quenched(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10_quenched_duration': 8.,'log10_total_duration': 9., 'log10Z': -3., 'log10M*': 0.})
sed = sed_modeller.get_Lnu(sfzh, {'fesc': 0.0}, dust = False)
ax.plot(sed.total.lam, np.log10(sed.total.lnu), c='0.5', lw=1)



# --- Binggeli
ax.axvline(3500., alpha=0.2, color='r', ls= '--', lw = 1, label = r'$\rm Binggeli+ $')
ax.axvline(4200., alpha=0.2, color='r', ls= '--', lw = 1)


# --- This work
ax.fill_between([3400, 3600], [0,0], [100,100], alpha=0.1, color='r', label = r'$\rm this\ work $')
ax.fill_between([4150, 4250], [0,0], [100,100], alpha=0.1, color='r')

# --- D4000
ax.fill_between([3750, 3950], [0,0], [100,100], alpha=0.1, color='g', label = r'$\rm D4000 $')
ax.fill_between([4050, 4250], [0,0], [100,100], alpha=0.1, color='g')




ax.legend(fontsize=7, labelspacing = 0.1)

ax.set_xlim([2500, 5000])

ax.set_ylim([18, 22])

ax.set_xlabel(r'$\rm \lambda/\AA$')
ax.set_ylabel(r'$\rm specific\ L_{\nu}/erg\ s^{-1}\ Hz^{-1}\ M_{\odot}^{-1}$')

fig.savefig(f'figs/sed_rest.pdf')
fig.clf()
