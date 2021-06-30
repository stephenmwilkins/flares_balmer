
import numpy as np
import matplotlib.cm as cm
import pickle


import flares
import flares_analysis
import flare.plt as fplt


import flare.filters

import h5py

# ----------------------------------------------------------------------
# --- open data


filedir = '/cosma7/data/dp004/dc-payy1/my_files/flares_pipeline/data/'

fl = flares.flares(filedir + 'flares.hdf5', sim_type='FLARES')
f = h5py.File(filedir + 'flares.hdf5', 'r')
# fl.explore()

model_tag = f'{fl.halos[0]}/{fl.tags[0]}/Galaxy/BPASS_2.2.1/Chabrier300'
lam = np.array(f[f'{model_tag}/SED/Wavelength'][:])



# for observatory in ['Webb','HubbleSpitzer']:
for observatory in ['HubbleSpitzer']:

    if observatory == 'Webb':
        filters = [f'Webb.NIRCAM.{filter}' for filter in ['F150W', 'F200W','F277W','F356W','F444W']]
    if observatory == 'HubbleSpitzer':
        filters = [f'Hubble.WFC3.{filter}' for filter in ['f105w','f125w','f160w']] + [f'Spitzer.IRAC.{filter}' for filter in ['ch1','ch2']]


    for log10Mstar_limit in [8.5, 9.5, 10.5]:
    # for log10Mstar_limit in [10.5]:

        # loop over simulations
        # for sim in fl.halos:

        O = {}

        # for i,tag in enumerate(fl.tags[0:3]):
        for i,tag in enumerate(fl.tags):

            redshift = fl.zeds[i]
            redshifts = np.arange(redshift-0.5, redshift+0.6, 0.1)

            O[tag] = {}
            O[tag]['z'] = redshifts

            for f1, f2 in zip(filters[:-1], filters[1:]):
                for p in [2.5, 16, 50, 84, 97.5]:
                    O[tag][f'{f1}_{f2}_P{p}'] = []


            for z in redshifts:

                lamz = lam*(1+z)

                # --- read in filter transmission functions and map to SED wavelength grid
                F = flare.filters.add_filters(filters, new_lam = lamz)

                O2 = {f:[] for f in filters}

                # for sim in [fl.halos[0]]:
                for sim in fl.halos:

                    model_tag = f'{sim}/{tag}/Galaxy/BPASS_2.2.1/Chabrier300'

                    log10Mstar_30 = np.log10(f[f'{sim}/{tag}/Galaxy/Mstar_30'][:])+10

                    s = log10Mstar_30>log10Mstar_limit

                    for spec_type in ['DustModelI']:
                        seds = np.array(f[f'{model_tag}/SED/{spec_type}'][s,:])

                    if len(seds)>0:
                        for sed in seds:
                            sed[lam<1216] = 0.0 # IGM
                            for filter in filters:
                                O2[filter].append(np.trapz(sed * F[filter].T, lamz) / np.trapz(F[filter].T, lamz))


                if len(O2[filters[0]])>0:
                    O[tag]['N'] = len(O2[filters[0]])
                    for f1, f2 in zip(filters[:-1], filters[1:]):
                        for P in [2.5, 16, 50, 84, 97.5]:
                            v = np.percentile(np.array(O2[f2])/np.array(O2[f1]), P)
                            O[tag][f'{f1}_{f2}_P{P}'].append(v)

                    print(f"{z:.1f} {len(O2[filters[0]])} {O[tag][f'{f1}_{f2}_P50'][-1]:.2f}")


        pickle.dump(O, open(f'data/{observatory}_colours_{log10Mstar_limit}.p', 'wb'))
