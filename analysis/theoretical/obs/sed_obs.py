




import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm

import interrogator
from interrogator.sed import models
import interrogator.sed.sfzh

import flare
import flare.plt as fplt
import flare.filters

cosmo = flare.default_cosmo()
z = 9.11
z = 10.

cf = lambda i: cm.rainbow(i/(len(filters)-1))

filters = ['Webb.NIRCam.F356W','Webb.NIRCam.F444W','Spitzer.IRAC.ch1', 'Spitzer.IRAC.ch2']
filters = ['Webb.NIRCam.F356W','Webb.NIRCam.F444W']
filters = ['Webb.NIRCam.F115W', 'Webb.NIRCam.F150W', 'Webb.NIRCam.F200W', 'Webb.NIRCam.F277W', 'Webb.NIRCam.F356W','Webb.NIRCam.F444W']


# -------------------------------------------------
# --- define choise of SPS model and initial mass function (IMF)

# SPS = models.SPS('P2/ModSalpeter_100')
sed_modeller = models.SPS('BPASSv2.2.1.binary/ModSalpeter_300')


# -------------------------------------------------
# --- define star formation and metal enrichment history (sfzh)



fig, ax = fplt.simple()

sfzh_constant, sfr = interrogator.sed.sfzh.constant(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10_duration': 8., 'log10Z': -3., 'log10M*': 8.})



sfzh_quenched, sfr = interrogator.sed.sfzh.quenched(sed_modeller.grid['log10age'], sed_modeller.grid['log10Z'] , {'log10_quenched_duration': 8.,'log10_total_duration': 9., 'log10Z': -3., 'log10M*': 8.})



#
# # --- D4000
# ax.fill_between([3750, 3950], [0,0], [100,100], alpha=0.1, color='k')
# ax.fill_between([4050, 4250], [0,0], [100,100], alpha=0.1, color='k')
#
#
# # --- Binggeli
# ax.axvline(3500., alpha=0.2, color='k', ls= '--', lw = 1)
# ax.axvline(4200., alpha=0.2, color='k', ls= '--', lw = 1)
#
#
# --- This work



BBa = np.array([3400, 3600])
BBb = np.array([4150, 4250])

ax.fill_between((1+z)*BBa/1E4, [-100,-100], [100,100], alpha=0.05, color='r')
ax.fill_between((1+z)*BBb/1E4, [-100,-100], [100,100], alpha=0.05, color='r')



#
# # --- get FUV and NUV luminosities
#
# filters = ['FAKE.FAKE.1500','FAKE.FAKE.2500']
# F = flare.filters.add_filters(filters, new_lam = SED.total.lam)
#
# L = SED.total.return_Lnu(F)
#
# print(L)


for sfzh in [sfzh_constant, sfzh_quenched]:

    sed = sed_modeller.get_Lnu(sfzh, {'fesc': 0.0}, dust = False)

    sed.total.get_fnu(cosmo, z) # calculate observer frame wavelength

    ax.plot(sed.total.lamz/1E4, np.log10(sed.total.fnu), zorder = 1, c='0.7', lw=1) # plot SED


    EW = sed.total.fnu[1:-1]/(0.5*(sed.total.fnu[0:-2]+sed.total.fnu[2:]))

    # ax.plot(sed.total.lamz[1:-1]/1E4, np.log10(EW))

    # --- this colour codes lines by their EW
    # norm =  mpl.colors.Normalize(vmin=1, vmax=3)
    #
    # for i,l in enumerate(sed.total.lamz[1:-1]):
    #     if EW[i]>10:
    #         l /= 1E4
    #         mn = np.log10(0.5*(sed.total.fnu[1:-1][i-1]+sed.total.fnu[1:-1][i+1]))
    #         mx = np.log10(sed.total.fnu[1:-1][i])
    #         ax.plot([l,l], [mn, mx], c=cm.magma(norm(np.log10(EW[i]))), lw=1)



    # --- add broadband photometry
    F = flare.filters.add_filters(filters, new_lam = sed.total.lamz) # --- NOTE: need to give it the redshifted
    sed.total.get_Fnu(F) # generates Fnu (broad band fluxes)
    for i,f in enumerate(filters):
        # ax.scatter(F[f].pivwv()/1E4, np.log10(sed.total.Fnu[f]), edgecolor = 'k', zorder = 2, label = f)
        ax.scatter(F[f].pivwv()/1E4, np.log10(sed.total.Fnu[f]), c = [cf(i)], zorder = 3, label = f, s=20)





for i,f in enumerate(filters):
    ax.plot(F[f].lam/1E4, F[f].T-1.01, c=cf(i))



ax.set_xlim([0.8, 5.5])

ax.set_ylim([-1,2.99])

ax.set_xlabel(r'$\rm \lambda_{obs}/\mu m$')
ax.set_ylabel(r'$\rm specific\ f_{\nu}/nJy\ (10^{8}\ M_{\odot})^{-1}$')

fig.savefig(f'figs/sed_obs.pdf')
fig.clf()
