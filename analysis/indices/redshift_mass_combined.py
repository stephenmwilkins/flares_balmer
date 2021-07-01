
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
log10Mstar_limit = 8.5


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


D = {}
s = {}
for tag, z in zip(fl.tags, fl.zeds):
    # --- get quantities (and weights and deltas)
    D[z] = fa.get_datasets(fl, tag, quantities)
    s[z] = D[z]['log10Mstar_30']>log10Mstar_limit


# --- get default limits and modify them to match the selection range
limits = fa.limits
limits['log10Mstar_30'] = [log10Mstar_limit, 10.9]

fig = fa.linear_redshift_mcol(D, fl.zeds, 'log10Mstar_30', properties, s, limits = limits, scatter_colour_quantity = 'log10FUV', scatter_cmap = cm.inferno)

fig.savefig(f'figs/combined_mass_redshift.pdf')
fig.savefig(f'figs/combined_mass_redshift.png')
