
import numpy as np
import matplotlib.cm as cm

import flares
import flares_analysis
import flare.plt as fplt

# ----------------------------------------------------------------------
# --- open data

fl = flares.flares('/cosma7/data/dp004/dc-payy1/my_files/flares_pipeline/data/flares.hdf5', sim_type='FLARES')

# fl.explore()

halo = fl.halos

# ----------------------------------------------------------------------
# --- define parameters and tag
tag = fl.tags[-3]  # --- select tag -3 = z=7
log10Mstar_limit = 8.5

q = 'beta'

# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []
quantities.append({'path': 'Galaxy', 'dataset': 'Mstar_30', 'name': None, 'log10': True})

for spec_type in ['Intrinsic','Pure_Stellar']:

    quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/beta', 'dataset': spec_type, 'name': f'beta_{spec_type}', 'log10': False})
    quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/BB_Wilkins', 'dataset': spec_type, 'name': f'BB_{spec_type}', 'log10': True})

    for f in ['FUV']:
        quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/{spec_type}', 'dataset': f, 'name': f'{f}_{spec_type}', 'log10': True})


# --- get quantities (and weights and deltas)
D = flares_analysis.get_datasets(fl, tag, quantities)



for q in ['log10BB','beta']:

    D['R'] = D[f'{q}_Intrinsic'] - D[f'{q}_Pure_Stellar']

    D[q] = D[f'{q}_Intrinsic']

    # ----------------------------------------------
    # define selection
    s = D['log10Mstar_30']>log10Mstar_limit


    # --- get default limits and modify them to match the selection range
    limits = flares_analysis.limits
    limits['log10Mstar_30'] = [log10Mstar_limit, 10.9]
    limits[q] = [np.percentile(D[q][s][D[q][s]==D[q][s]], 1), np.max(D[q][s][D[q][s]==D[q][s]])]
    limits['R'] = [np.percentile(D['R'][s][D['R'][s]==D['R'][s]], 1), np.max(D['R'][s][D['R'][s]==D['R'][s]])]

    labels = flares_analysis.labels
    labels['R'] = rf'{labels[q]}_{{intrinsic}} - {labels[q]}_{{pure\ stellar}}'
    labels[q] = rf'{labels[q]}_{{intrinsic}}'


    # --- this is useful
    fig, ax, cax = flares_analysis.simple_wcbar(D, 'log10Mstar_30', 'R', q, s, limits = limits, labels = labels, cmap = cm.cividis, add_weighted_median= True)

    fig.savefig(f'figs/{q}_nebular.pdf')
    fig.clf()
