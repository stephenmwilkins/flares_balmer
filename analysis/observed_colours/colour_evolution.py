
import matplotlib as mpl
import matplotlib.cm as cm
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import numpy as np

import pickle

import flares
import flares_analysis
import flare.plt as fplt




observatory = 'HubbleSpitzer'
observatory = 'Webb'




for observatory in ['Webb','HubbleSpitzer']:

    if observatory == 'Webb':
        filters = [f'Webb.NIRCAM.{filter}' for filter in ['F150W', 'F200W','F277W','F356W','F444W']]
    if observatory == 'HubbleSpitzer':
        filters = [f'Hubble.WFC3.{filter}' for filter in ['f105w','f125w','f160w']] + [f'Spitzer.IRAC.{filter}' for filter in ['ch1','ch2']]


    O = {}
    for log10Mstar_limit in [8.5, 9.5, 10.5]:
        O[log10Mstar_limit] = pickle.load(open(f'data/{observatory}_colours_{log10Mstar_limit}.p','rb'))


    print(list(O.keys()))


    # --- define colour scale
    norm = mpl.colors.Normalize(vmin=5, vmax=10)

    fig, axes = plt.subplots(len(filters)-1, 1, figsize = (3.5,5), sharex = True)
    plt.subplots_adjust(left=0.15, top=0.95, bottom=0.15, right=0.9, wspace=0.0, hspace=0.0)


    for f1, f2, ax in zip(filters[:-1], filters[1:], axes):

        for tag, z in zip(O[8.5].keys(), np.arange(10, 4, -1)):

            c = cm.plasma(norm(z))

            ax.fill_between(O[8.5][tag]['z'], np.log10(O[8.5][tag][f'{f1}_{f2}_P16']), np.log10(O[8.5][tag][f'{f1}_{f2}_P84']), color=c, alpha=0.2)
            ax.plot(O[8.5][tag]['z'], np.log10(O[8.5][tag][f'{f1}_{f2}_P50']), c=c, lw=1)

            ax.plot(O[9.5][tag]['z'], np.log10(O[9.5][tag][f'{f1}_{f2}_P50']), c=c, lw=1, ls='-.')
            if len(O[10.5][tag][f'{f1}_{f2}_P50']) == len(O[10.5][tag]['z']):
                ax.plot(O[10.5][tag]['z'], np.log10(O[10.5][tag][f'{f1}_{f2}_P50']), c=c, lw=1, ls=':')

        ax.set_xlim(4.5, 10.5)
        ax.set_ylim(-0.3, 0.5)
        ax.set_ylabel(rf"$\rm log_{{10}}(f_{{ {f2.split('.')[-1]} }}/f_{{ {f1.split('.')[-1]} }}) $", fontsize=8)


        if f1=='Spitzer.IRAC.ch1' and f2=='Spitzer.IRAC.ch2':

            observations = []
            observations.append({'id': 'MACS1149-JD1', 'ref': 'Hashimoto', 'z': 9.1096, 'F444W-F356W': [0.92,0.19]})
            observations.append({'id': 'GN-z9-1', 'ref': 'Laporte+2021', 'z': 9.89, 'F444W-F356W': [0.24, 0.33]})
            observations.append({'id': 'GN-z10-3', 'ref': 'Laporte+2021', 'z': 8.78, 'F444W-F356W': [0.24, 0.33]})
            # observations.append({'id': 'GS-z9-1', 'ref': 'Laporte+2021', 'z': , 'F444W-F356W': [0.89, 0.49]})
            observations.append({'id': 'MACS0416-JD', 'ref': 'Laporte+2021', 'z': 9.28, 'F444W-F356W': [0.54, None]})
            # observations.append({'id': 'UVISTA-1212', 'ref': 'Laporte+2021', 'z': , 'F444W-F356W': [1.05, 0.47]})


            for obs, ms in zip(observations, ['o','d','h','p']):

                m = 0.4*obs['F444W-F356W'][0]

                axes[-1].scatter(obs['z'], m , c='k', zorder = 3, marker = ms, label = rf'$\rm {obs["id"]}$', s=15)

                if obs['F444W-F356W'][1]:
                    e = [m-0.4*obs['F444W-F356W'][1], m+0.4*obs['F444W-F356W'][1]]
                    axes[-1].plot([obs['z']]*2, e , c='k', lw=1)
                else:
                    # axes[-1].arrow(, 0.0, 0.05, head_width = 0.1, color='k')
                    axes[-1].annotate("", xy=(obs['z'], m+0.2), xytext=(obs['z'], m), arrowprops=dict(arrowstyle="->"))


    axes[-1].set_xlabel(rf'$\rm z $')


    handles = [Line2D([0], [0], label = r'$\rm log_{10}(M_{\star}/M_{\odot})>8.5$' , color='k', lw=1, ls='-')]
    handles += [Line2D([0], [0], label = r'$\rm log_{10}(M_{\star}/M_{\odot})>9.5$' , color='k', lw=1, ls='-.')]
    handles += [Line2D([0], [0], label = r'$\rm log_{10}(M_{\star}/M_{\odot})>10.5$' , color='k', lw=1, ls=':')]


    axes[0].legend(handles=handles, fontsize=7, labelspacing = 0.1)
    axes[-1].legend(fontsize=5, labelspacing = 0.1)

    fig.savefig(f'figs/colour_evolution_{observatory}.pdf')
