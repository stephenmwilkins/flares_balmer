
import numpy as np
import matplotlib.cm as cm

import flares
import flares_analysis
import flare.plt as fplt

# ----------------------------------------------------------------------
# --- open data

fl = flares.flares('/cosma7/data/dp004/dc-payy1/my_files/flares_pipeline/data/flares.hdf5', sim_type='FLARES')

# fl.explore()

# ----------------------------------------------------------------------
# --- define parameters and tag
tag = fl.tags[-3]  # --- select tag -3 = z=7
log10Mstar_limit = 8.5


# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []
quantities.append({'path': 'Galaxy', 'dataset': 'SFR_inst_30', 'name': None, 'log10': True})
quantities.append({'path': 'Galaxy', 'dataset': 'Mstar_30', 'name': None, 'log10': True})
# need age


quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/beta', 'dataset': 'Intrinsic', 'name': 'beta', 'log10': False})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/BB_Wilkins', 'dataset': 'Intrinsic', 'name': 'BB', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/Intrinsic/HI4861', 'dataset': 'EW', 'name': 'HbetaEW', 'log10': True})


for spec_type in ['Intrinsic','DustModelI']:
    for f in ['FUV', 'V']:
        quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/{spec_type}', 'dataset': f, 'name': f'{f}_{spec_type}', 'log10': True})





# --- get quantities (and weights and deltas)
D = flares_analysis.get_datasets(fl, tag, quantities)


# ----------------------------------------------
# define new quantities
D['log10sSFR'] = D['log10SFR_inst_30']-D['log10Mstar_30']+9
D['AFUV'] = 2.5*(D['log10FUV_Intrinsic']-D['log10FUV_DustModelI'])
D['log10M*/LV'] = D['log10Mstar_30']-D['log10V_DustModelI']


# ----------------------------------------------
# define selection
s = D['log10Mstar_30']>log10Mstar_limit

# ----------------------------------------------
# Print number of galaxies meeting the selection
print(f"Total number of galaxies: {len(D['log10Mstar_30'][s])}")


for q in ['log10BB', 'log10HbetaEW','beta']:

    # --- define quantities for including in the corner plot
    properties = [q, 'log10sSFR', 'log10sSFR']

    # --- get default limits and modify them to match the selection range
    limits = flares_analysis.limits
    limits['log10Mstar_30'] = [log10Mstar_limit, 10.9]

    labels = flares_analysis.labels
    labels[q] = f'{labels[q]}_{{intrinsic}}'



    # --- make corner plot
    fig = flares_analysis.linear(D, properties, s, limits = limits, scatter_colour_quantity = 'log10Mstar_30', scatter_cmap = cm.viridis, full_width = False)

    fig.savefig(f'figs/{q}_intrinsic_properties.pdf')
