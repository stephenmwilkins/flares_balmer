
import numpy as np
import matplotlib.cm as cm

import flares
import flares_analysis as fa
import flare.plt as fplt

# ----------------------------------------------------------------------
# --- open data

fl = flares.flares('/cosma7/data/dp004/dc-payy1/my_files/flares_pipeline/data/flares.hdf5', sim_type='FLARES')

# fl.explore()

halo = fl.halos

q = 'beta'


# ----------------------------------------------------------------------
# --- define parameters and tag
tag = fl.tags[-3]  # --- select tag -3 = z=7


s_limit = {'log10Mstar_30': 8.5, 'log10FUV': 28.5}

# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []

quantities.append({'path': 'Galaxy', 'dataset': 'Mstar_30', 'name': None, 'log10': True})
quantities.append({'path': 'Galaxy', 'dataset': 'SFR_inst_30', 'name': None, 'log10': True})

# for spec_type in ['Intrinsic','DustModelI']:
#     for f in ['FUV']:
#         quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/{spec_type}', 'dataset': f, 'name': f'{f}_{spec_type}', 'log10': True})

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI', 'dataset': 'FUV', 'name': None, 'log10': True})

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/beta', 'dataset': 'DustModelI', 'name': 'beta', 'log10': False})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/BB_Wilkins', 'dataset': 'DustModelI', 'name': 'BB', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/HI4861', 'dataset': 'EW', 'name': 'HbetaEW', 'log10': True})




properties = ['log10BB', 'log10HbetaEW','beta']
x_ = ['log10FUV', 'log10Mstar_30']

D = {}
s = {}
s['log10Mstar_30'] = {}
s['log10FUV'] = {}
for tag, z in zip(fl.tags, fl.zeds):

    # --- get quantities (and weights and deltas)
    D[z] = fa.get_datasets(fl, tag, quantities)

    for x in x_:
        s[x][z] = D[z][x]>s_limit[x]



# --- get default limits and modify them to match the selection range

cmap = {'log10FUV': cm.viridis, 'log10Mstar_30': cm.inferno}




for x, z in zip(x_, x_[::-1]): #the colour map is for the other parameter

    limits = fa.limits
    limits[x][0] = s_limit[x]

    fig = fa.linear_redshift_mcol(D, fl.zeds, x, properties, s[x], limits = limits, scatter_colour_quantity = z, scatter_cmap = cmap[x], add_linear_fit = True)

    fig.savefig(f'figs/combined_redshift_{x}.pdf')
    fig.savefig(f'figs/combined_redshift_{x}.png')
