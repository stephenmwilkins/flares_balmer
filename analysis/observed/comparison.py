
import numpy as np
import matplotlib.cm as cm

import flares
import flares_analysis
import flares_utility.analyse as analyse
import flare.plt as fplt
import flares_utility.limits
import flares_utility.plt
import flares_utility.colors
# ----------------------------------------------------------------------
# --- open data

filename = '/Users/stephenwilkins/Dropbox/Research/data/simulations/flares/flares_noparticlesed.hdf5'
flares = analyse.analyse(filename, default_tags = False)



# ----------------------------------------------------------------------
# --- define parameters and tag
tag = flares.tags[-3]  # --- select tag -3 = z=7
z = 'log10FUV'

# ----------------------------------------------------------------------
# --- define quantities to read in [not those for the corner plot, that's done later]

quantities = []

quantities.append({'path': 'Galaxy/Mstar_aperture', 'dataset': f'30', 'name': 'Mstar_30', 'log10': True})

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/beta', 'dataset': f'DustModelI', 'name': f'beta', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/BB_Wilkins', 'dataset': 'DustModelI', 'name': 'BB', 'log10': True})

# quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/HI4861', 'dataset': 'EW', 'name': 'HbetaEW', 'log10': True})

quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/HI4861', 'dataset': 'EW', 'name': 'HI4861_EW', 'log10': True})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/OIII5007', 'dataset': 'EW', 'name': 'OIII5007_EW', 'log10': False})
quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Lines/DustModelI/OIII4959', 'dataset': 'EW', 'name': 'OIII4959_EW', 'log10': False})


quantities.append({'path': f'Galaxy/BPASS_2.2.1/Chabrier300/Luminosity/DustModelI', 'dataset': 'FUV', 'name': None, 'log10': True})




# --- get quantities (and weights and deltas)
D = flares.get_datasets(tag, quantities)

D['OIIIHbEW'] = D['HI4861_EW']+D['OIII5007_EW']+D['OIII4959_EW']
D['log10OIIIHbEW'] = np.log10(D['OIIIHbEW'])



# ----------------------------------------------
# define selection
s = D['log10FUV']>28.

# ----------------------------------------------
# Print number of galaxies meeting the selection
print(f"Total number of galaxies: {np.sum(s)}")


# ----------------------------------------------
# ----------------------------------------------
# Corner Plot

# --- define quantities for including in the corner plot
properties = ['beta', 'log10BB', 'log10OIIIHbEW']

# --- get default limits and modify them to match the selection range
limits = flares_utility.limits.limits



# --- make corner plot
fig = flares_utility.plt.corner(D, properties, s, limits = limits, scatter_colour_quantity = z, scatter_cmap = flares_utility.colors.cmap[z], full_width = False)

# fig = flares_analysis.corner3(D, properties, s, {'beta': cm.viridis, 'log10BB': cm.plasma, 'log10HbetaEW': cm.inferno}, limits = limits)


fig.savefig(f'figs/comparison.pdf')
