
import numpy as np

import matplotlib as mpl
import matplotlib.cm as cm
mpl.use('Agg')
import matplotlib.pyplot as plt

import flares
import flares_analysis as fa
import flare.plt as fplt

# ----------------------------------------------------------------------
# --- open data

fl = flares.flares('/cosma7/data/dp004/dc-payy1/my_files/flares_pipeline/data/flares.hdf5', sim_type='FLARES')


# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []

quantities.append({'path': 'Galaxy', 'dataset': 'Mstar_30', 'name': None, 'log10': True})

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/beta', 'dataset': 'DustModelI', 'name': 'beta', 'log10': False})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/BB_Wilkins', 'dataset': 'DustModelI', 'name': 'BB', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/HI4861', 'dataset': 'EW', 'name': 'HbetaEW', 'log10': True})


for q, c in zip(['log10BB', 'log10HbetaEW','beta'], ['r','g','b']):

    O = {}
    O['z'] = fl.zeds

    log10Mstar_limits = [8.5, 9., 9.5, 10.]

    for log10Mstar_limit in log10Mstar_limits:

        O[log10Mstar_limit] = {}

        for p in [2.5, 16, 50, 84, 97.5]:
            O[log10Mstar_limit][p] = []

    for tag, z in zip(fl.tags, fl.zeds):

        # --- get quantities (and weights and deltas)
        D = fa.get_datasets(fl, tag, quantities)

        for log10Mstar_limit in log10Mstar_limits:

            # s = (D['log10Mstar_30']>log10Mstar-0.25)&(D['log10Mstar_30']<log10Mstar+0.25)
            s = D['log10Mstar_30']>log10Mstar_limit

            for p in [2.5, 16, 50, 84, 97.5]:
                O[log10Mstar_limit][p].append(np.percentile(D[q][s], p))



    fig = plt.figure(figsize = (3.5, 2.5))

    left  = 0.15
    height = 0.8
    bottom = 0.15
    width = 0.8

    ax = fig.add_axes((left, bottom, width, height))

    ax.fill_between(fl.zeds, O[8.5][2.5], O[8.5][97.5], alpha= 0.05, color = c)
    ax.fill_between(fl.zeds, O[8.5][16], O[8.5][84], alpha= 0.05, color = c)

    for log10Mstar_limit, ls in zip(log10Mstar_limits, ['-','--','-.',':']):
        ax.plot(fl.zeds, O[log10Mstar_limit][50], ls = ls, c=c, lw=1, label = r'$\rm log_{10}(M_{\star}/M_{\odot})>'+str(log10Mstar_limit)+'$')

    ax.set_ylabel(rf'$\rm {fa.labels[q]}$', fontsize = 9)
    ax.set_xlabel(rf'$\rm z$', fontsize = 9)

    ax.legend(fontsize = 8)


    fig.savefig(f'figs/{q}_redshift_summary.pdf')
