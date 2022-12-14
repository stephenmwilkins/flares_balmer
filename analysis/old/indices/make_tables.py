
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
log10FUV_limit = 28.5

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
s['log10Mstar_30'] = {}
s['log10FUV'] = {}


for tag, z in zip(fl.tags, fl.zeds):
    # --- get quantities (and weights and deltas)
    D[z] = fa.get_datasets(fl, tag, quantities)
    s['log10Mstar_30'][z] = D[z]['log10Mstar_30']>log10Mstar_limit
    s['log10FUV'][z] = D[z]['log10FUV']>log10FUV_limit




norm = {'log10Mstar_30': 9.,'log10FUV': 29.}


for j,y in enumerate(properties):

    print(r'\hline')
    print(rf'\multicolumn{{9}}{{c}}{{$\rm {fa.labels[y]} $}} \\')
    print(r'\hline')

    for i,z in enumerate(fl.zeds):

        fit_p = {}
        rho = {}
        N = {}


        for x in ['log10Mstar_30','log10FUV']:

            N[x] = len(D[z]['weight'][s[x][z]])

            x_ = D[z][x][s[x][z]] - norm[x]
            y_ = D[z][y][s[x][z]]
            w_ = D[z]['weight'][s[x][z]]
            s_ = (~np.isnan(x_))&(~np.isnan(y_))&(~np.isinf(y_))

            fit_p[x] = np.polyfit(x_[s_], y_[s_], 1, w = w_[s_])

            rho[x] = np.corrcoef(x_[s_], y_[s_])



        # if i == 0:
        #     print(rf" \multirow{{ {len(fl.zeds)} }}{{4em}}{{ {fa.labels[y]} }}& {z} & {fit_p['log10Mstar_30'][1]:.2f} & {fit_p['log10Mstar_30'][0]:.2f} & {fit_p['log10FUV'][1]:.2f} & {fit_p['log10FUV'][0]:.2f} \\")
        # else:
        #     print(rf" & {z} & {fit_p['log10Mstar_30'][1]:.2f} & {fit_p['log10Mstar_30'][0]:.2f} & {fit_p['log10FUV'][1]:.2f} & {fit_p['log10FUV'][0]:.2f} \\")

        print(rf" {int(z)} & {N['log10Mstar_30']} & {fit_p['log10Mstar_30'][1]:.2f} & {fit_p['log10Mstar_30'][0]:.2f} & {rho['log10Mstar_30'][0][1]:.2f} & {N['log10FUV']} & {fit_p['log10FUV'][1]:.2f} & {fit_p['log10FUV'][0]:.2f} & {rho['log10FUV'][0][1]:.2f} \\")




# \multirow{3}{4em}{Multiple row} & cell2 & cell3 \\
# & cell5 & cell6 \\
# & cell8 & cell9 \\
# \hline
