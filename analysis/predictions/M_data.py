
import numpy as np
from astropy.table import Table, Column

import flares_utility.analyse as analyse
import flares_utility.stats as stats

from unyt import unyt_quantity

tags = ['005_z010p000', '006_z009p000', '007_z008p000',
        '008_z007p000', '009_z006p000', '010_z005p000']
redshifts = [10, 9, 8, 7, 6, 5]


filename = '/Users/sw376/Dropbox/Research/data/simulations/flares/flares_noparticlesed_v4.hdf5'

flares = analyse.analyse(filename, default_tags=False)

flares.list_datasets()

binw = 0.2
quantiles = [0.022, 0.158, 0.50, 0.842, 0.978]
percentiles = [np.round(q*100, 1) for q in quantiles]



qname = 'BB'
qpath = f'Galaxy/BPASS_2.2.1/Chabrier300/Indices/BB_Wilkins'
qds = 'DustModelI'

unit_str = 'dimensionless'
description = False
logged = True


print('-'*5, qname, unit_str, description)

# ---
quantities = []
quantities.append({'path': 'Galaxy/Mstar_aperture',
                    'dataset': f'30', 'name': 'Mstar', 'log10': True})
quantities.append({'path': qpath,
                    'dataset': qds, 'name': 'Y', 'log10': logged})

t = Table()

z_ = np.array([])
bin_centres_ = np.array([])
N_ = np.array([])
y_ = {}

for q in quantiles:
    y_[q] = np.array([])

for z in redshifts:

    tag = flares.tag_from_zed[z]

    D = flares.get_datasets(tag, quantities)
    X = D['log10Mstar']
    w = D['weight']

    if logged:
        Y = D['log10Y']
    else:
        Y = D['Y']

    bins = np.arange(8., np.max(X), binw)
    bin_centres = 0.5*(bins[:-1]+bins[1:])

    out = stats.binned_weighted_quantile(X, Y, w, bins, quantiles)
    N, _ = np.histogram(X, bins=bins)

    # --- make stacks
    z_ = np.hstack((z_, np.ones(len(N))*z))
    bin_centres_ = np.hstack((bin_centres_, bin_centres))
    N_ = np.hstack((N_, N))
    for i, q in enumerate(quantiles):
        y_[q] = np.hstack((y_[q], np.array(out[:, i])))

t.add_column(Column(data=z_, name='z'))
t.add_column(Column(data=np.round(bin_centres_, 2), name='log10Mstar', unit='dex(solMass)',
                description=r'log10(stellar mass)'))
t.add_column(Column(data=N_.astype(int), name='N'))

for q in quantiles:
    if logged:
        t.add_column(Column(data=np.round(
            y_[q], 2), name=f'log10{qname}_P{q*100:.1f}', unit='dex('+unit_str+')', description=description))
    else:
        t.add_column(Column(data=np.round(
            y_[q], 2), name=f'{qname}_P{q*100:.1f}', unit=unit_str, description='log10('+description+')'))

t.meta['name'] = 'FLARES'
t.meta['references'] = []
t.meta['percentiles'] = percentiles

t.write(
    f'data/Mstar.ecsv', format='ascii.ecsv', overwrite=True)
