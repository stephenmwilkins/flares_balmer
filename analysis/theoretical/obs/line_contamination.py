

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import flare
import flare.plt as fplt
import flare.filters

import interrogator
from interrogator.sed import models
import interrogator.sed.sfzh




# -------------------------------------------------
# --- define choise of SPS model and initial mass function (IMF)


line_modeller = models.lines('BPASSv2.2.1.binary/ModSalpeter_300')

# -------------------------------------------------
# --- define star formation and metal enrichment history (sfzh)

sfzh, sfr = interrogator.sed.sfzh.constant(line_modeller.grid['log10age'], line_modeller.grid['log10Z'] , {'log10_duration': 7., 'log10Z': -2., 'log10M*': 8.})

fesc = 0.0
li = line_modeller.get_all_info_cols(sfzh, {'fesc': fesc}, dust_model = ('just_gas'))

print(li)

filters = ['Webb.NIRCam.F115W', 'Webb.NIRCam.F150W', 'Webb.NIRCam.F200W', 'Webb.NIRCam.F277W', 'Webb.NIRCam.F356W', 'Webb.NIRCam.F410M', 'Webb.NIRCam.F430M', 'Webb.NIRCam.F444W']
filters = ['Webb.NIRCam.F115W', 'Webb.NIRCam.F150W', 'Webb.NIRCam.F200W', 'Webb.NIRCam.F277W', 'Webb.NIRCam.F356W','Webb.NIRCam.F444W']
# filters = ['Webb.NIRCam.F444W']
F = flare.filters.add_filters(filters)


s = (li['lam']*(1+12)>18000.)&(li['lam']*(1+5)<50000.)&(li['EW']>10.)




fig, axes = plt.subplots(6, 1, sharex = True, figsize = (3.5, 2.0))
plt.subplots_adjust(bottom=0.15, top=0.95, wspace=0.0, hspace=0.0)


z = np.arange(5,12,0.01)

for iax, (f, ax) in enumerate(zip(filters, axes)):
    print('-'*20)
    print(f, F[f].mnmx())
    mn, mx = F[f].mnmx()
    piv = F[f].pivwv()

    line_EW = np.zeros(z.shape)



    # --- determine line contamination
    for i, (line, lam, EW) in enumerate(zip(li['line'][s],li['lam'][s], li['EW'][s])):
        mnz, mxz = mn/lam-1., mx/lam-1.
        c = cm.magma(i/len(li['line'][s]))
        line_EW[(z<mxz)&(z>mnz)] += EW

        if mnz < 12. and mxz > 5.:
            ax.fill_between([mnz, mxz], [0,0], [1,1], alpha = EW/2000., color=c)



    # --- determine UVC suitability

    mnz, mxz = mn/1216-1, mx/3700-1 # defines that a filter is entirely above the Balmer break
    # ax.fill_between([mnz, mxz], [0,0], [1,1], facecolor="none", edgecolor='0.8', hatch='...')
    ax.plot([mnz,mxz],[0.5,0.5], c='blue', lw=5, alpha=0.3)

    # --- determine BB suitability

    # mnz, mxz = mn/3700-1, mx/5000-1 # defines that a filter is entirely above the Balmer break
    # ax.fill_between([mnz, mxz], [0,0], [1,1], facecolor="none", edgecolor='0.5', hatch='//')

    mnz, mxz = piv/2500-1, piv/3700-1 # defines that a filter is entirely above the Balmer break
    # ax.fill_between([mnz, mxz], [0,0], [1,1], facecolor="none", edgecolor='0.5', hatch='\\\\\\\\')
    # ax.axvline(mnz)
    # ax.axvline(mxz)
    ax.plot([mnz,mxz],[0.5,0.5], c='b', alpha=1.0)

    mxz, mnz = mn/3700-1, mx/7000-1 # defines that a filter is entirely above the Balmer break
    mxz, mnz = piv/3700-1, piv/7000-1 # defines that a filter is entirely above the Balmer break
    #ax.fill_between([mnz, mxz], [0,0], [1,1], facecolor="none", edgecolor='0.5', hatch='////')
    # ax.axvline(mnz)
    # ax.axvline(mxz)

    ax.plot([mnz,mxz],[0.5,0.5], c='r', alpha=0.5)

    mnz2 = z[(z>mnz)&(line_EW<200)][0]

    print('good BB red', f, np.round(mnz2,1), np.round(mxz,1))

    ax.plot([mnz2,mxz],[0.5,0.5], c='r', alpha=1.0)

    if iax>0:
        ax.plot([mnz2,mnz2,mxz,mxz],[1,0.01,0.01,1], c='k', alpha=1.0, lw=1)
        axes[iax-1].plot([mnz2,mnz2,mxz,mxz],[0,0.99,0.99,0], c='k', alpha=1.0, lw=1)


    ax.text(4.9, 0.5, f.split('.')[-1], va='center', ha='right', fontsize=6)

    ax.set_ylim([0, 1])
    ax.set_yticks([])
    ax.set_xlim([5., 12])

    if ax == axes[-1]: ax.set_xlabel(r'$\rm z$')

fig.savefig(f'figs/line_contamination.pdf')
fig.clf()














# ----- second plot


fig = plt.figure(figsize = (3.5, 3))

left  = 0.15
height = 0.8
bottom = 0.15
width = 0.7

ax = fig.add_axes((left, bottom, width, height))


# --- determine line contamination
for i, (line, lam, EW) in enumerate(zip(li['line'][s],li['lam'][s], li['EW'][s])):
    # c = cm.magma(i/len(li['line'][s]))
    c = 'k'
    ax.axhline(lam, lw=1, c=c, alpha = EW/1000)
    ax.text(12.1, lam, rf"{line}", va='center', fontsize=5, c=c, alpha = EW/1000)

# --- plot break boundary

ax.axhline(1216, c='k', lw=1, ls='--', alpha=1.0)
ax.axhline(3700, c='k', lw=1, ls='--', alpha=1.0)

z = np.arange(5, 12, 0.01)

for i,f in enumerate(filters):
    mn, mx = F[f].mnmx()
    piv = F[f].pivwv()
    c = cm.rainbow(i/len(filters))

    d = mn/(1+z)
    u = mx/(1+z)

    ax.fill_between(z, d, u, color=c, alpha = 0.5)
    ax.text(12.1, piv/13., f.split('.')[-1], fontsize=5, c=c)


# ax.set_ylim([0, 1])
# ax.set_xlim([5., 12])
ax.set_xlim([5., 12])
ax.set_xlabel(r'$\rm z$')
ax.set_ylabel(r'$\rm \lambda_{rest}/\AA$')

fig.savefig(f'figs/line_contamination2.pdf')
fig.clf()
