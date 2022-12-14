
import numpy as np
import matplotlib.cm as cm

import flares
import flares_analysis
import flare.plt as fplt


# ----------------------------------------------------------------------
# --- open data

fl = flares.flares('/cosma7/data/dp004/dc-payy1/my_files/flares_pipeline/data/flares.hdf5', sim_type='FLARES')

halo = fl.halos

# ----------------------------------------------------------------------
# --- define parameters and tag
tag = fl.tags[-3]  # --- select tag -3 = z=7
log10Mstar_limit = 8.5



# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []
quantities.append({'path': 'Galaxy', 'dataset': 'Mstar_30', 'name': None, 'log10': True})

for spec_type in ['Intrinsic','DustModelI']:

    quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/beta', 'dataset': spec_type, 'name': f'beta_{spec_type}', 'log10': False})
    quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/BB_Wilkins', 'dataset': spec_type, 'name': f'BB_{spec_type}', 'log10': True})
    quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/{spec_type}/HI4861', 'dataset': 'EW', 'name': f'HbetaEW_{spec_type}', 'log10': True})

    for f in ['FUV']:
        quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/{spec_type}', 'dataset': f, 'name': f'{f}_{spec_type}', 'log10': True})


# --- get quantities (and weights and deltas)
D = flares_analysis.get_datasets(fl, tag, quantities)

# ----------------------------------------------
# define new quantities
D['AFUV'] = 2.5*(D['log10FUV_Intrinsic']-D['log10FUV_DustModelI'])



for q in ['log10BB', 'log10HbetaEW','beta']:

    D['R'] = D[f'{q}_DustModelI'] - D[f'{q}_Intrinsic']

    D[q] = D[f'{q}_DustModelI']

    # ----------------------------------------------
    # define selection
    s = D['log10Mstar_30']>log10Mstar_limit


    # --- get default limits and modify them to match the selection range
    limits = flares_analysis.limits
    limits['log10Mstar_30'] = [log10Mstar_limit, 10.9]
    limits['R'] = [np.percentile(D['R'][s][D['R'][s]==D['R'][s]], 1), np.max(D['R'][s][D['R'][s]==D['R'][s]])]

    labels = flares_analysis.labels

    labels['R'] = rf'{labels[q]} - {labels[q]}_{{intrinsic}}'

    # --- this is useful
    fig, ax, cax = flares_analysis.simple_wcbar(D, 'log10Mstar_30', 'R', q, s, limits = limits, labels = labels, cmap = cm.cividis, add_weighted_median= False)

    fig.savefig(f'figs/{q}_dust.pdf')
    fig.clf()


    # --- this is useful
    fig, ax, cax = flares_analysis.simple_wcbar(D, 'AFUV', 'R', 'log10Mstar_30', s, limits = limits, labels = labels, cmap = cm.viridis, add_weighted_median= False)

    fig.savefig(f'figs/{q}_dust1.pdf')
    fig.clf()
