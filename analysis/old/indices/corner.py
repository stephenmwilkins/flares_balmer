
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
quantities.append({'path': 'Galaxy', 'dataset': 'SFR_inst_30', 'name': None, 'log10': True})
quantities.append({'path': 'Galaxy', 'dataset': 'Mstar_30', 'name': None, 'log10': True})



for spec_type in ['Intrinsic','DustModelI']:

    quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/beta', 'dataset': spec_type, 'name': f'beta_{spec_type}', 'log10': False})

    for f in ['NUV', 'FUV']:
        quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/{spec_type}', 'dataset': f, 'name': f'{spec_type}_{f}', 'log10': True})





# --- get quantities (and weights and deltas)
D = flares_analysis.get_datasets(fl, tag, quantities)


# ----------------------------------------------
# define new quantities
D['log10sSFR'] = np.log10(D['SFR_inst_30'])-np.log10(D['Mstar_30'])+9

for spec_type in ['Intrinsic','DustModelI']:
    D[f'beta2_{spec_type}'] = np.log10(D[spec_type+'_FUV']/D[spec_type+'_NUV'])/np.log10(1500/2500)-2.0


# ----------------------------------------------
# define selection
s = D['log10Mstar_30']>log10Mstar_limit

# ----------------------------------------------
# Print number of galaxies meeting the selection
print(f"Total number of galaxies: {len(D['log10Mstar_30'][s])}")


# ----------------------------------------------
# ----------------------------------------------
# Corner Plot

# --- define quantities for including in the corner plot
properties = ['log10Mstar_30', 'log10sSFR', 'beta_DustModelI', 'beta2_DustModelI']

# --- get default limits and modify them to match the selection range
limits = flares_analysis.limits
limits['log10Mstar_30'] = [log10Mstar_limit, 10.9]


# --- make corner plot
fig = flares_analysis.corner_plot(D, properties, s, limits = limits, scatter_colour_quantity = 'log10Mstar_30', scatter_cmap = cm.viridis)



fig.savefig(f'figs/corner.pdf')
