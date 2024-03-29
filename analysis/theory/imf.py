

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.cm as cm
import cmasher as cmr
from astropy.table import Table

import sys
import os

from synthesizer.plt import single
from synthesizer.filters import TopHatFilterCollection
from synthesizer.grid import SpectralGrid, parse_grid_id
from synthesizer.parametric.sfzh import SFH, ZH, generate_sfzh, generate_instant_sfzh
from synthesizer.parametric.galaxy import SEDGenerator
from unyt import yr, Myr


if __name__ == '__main__':

    # -------------------------------------------------
    # --- calcualte the EW for a given line as a function of age

    models = [
        'bpass-v2.2.1-bin_100-300_cloudy-v17.03_log10Uref-2',
        'bpass-v2.2.1-bin_135-300_cloudy-v17.03_log10Uref-2',
        'bpass-v2.2.1-bin_135-100_cloudy-v17.03_log10Uref-2',
        'bpass-v2.2.1-bin_170-300_cloudy-v17.03_log10Uref-2',
    ]

    Z = 0.005
    log10Z = np.log10(Z)

    Zh = ZH.deltaConstant({'Z': Z})  # constant metallicity

    fig, ax = single()

    for model, ls in zip(models, ['-.', '-', ':', '--']):

        model_info = parse_grid_id(model)
        print(model_info)

        grid = SpectralGrid(model)

        log10durations = np.arange(0., 3.1, 0.1)

        label = model_info['imf'] + rf'$\ \rm m_{{\rm up}}={model_info["imf_hmc"]}$'

        log10bb = []
        age = []

        for log10duration in log10durations:

            # --- define the functional form of the star formation and metal enrichment histories
            sfh = SFH.Constant({'duration': 10**log10duration * Myr})

            # --- get the 2D star formation and metal enrichment history for the given SPS grid. This is (age, Z).
            sfzh = generate_sfzh(grid.log10ages, grid.metallicities, sfh, Zh)

            age_ = sfzh.calculate_median_age()/1E6
            age.append(age_)
            # --- define galaxy object
            # by default this automatically calculates the pure stellar spectra
            galaxy = SEDGenerator(grid, sfzh)
            galaxy.pacman()  # adds nebular emission

            # --- get quanitities

            sed = galaxy.spectra['total']
            log10bb_ = sed.get_balmer_break()
            log10bb.append(log10bb_)

        ax.plot(np.log10(age), log10bb, lw=1, color='k', alpha=1, ls=ls, label=label)

    ax.set_xlim([0, 3.])
    ax.set_ylabel(r'$\rm log_{10}(L_{4200}/L_{3500})$')
    ax.legend(fontsize=7, labelspacing=0.1)
    ax.set_xlabel(r'$\rm log_{10}(age/Myr)$')

    fig.savefig('figs/theory_imf.pdf')
