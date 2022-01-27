# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import os
import numpy as np

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_ent_N_flow(_model, _file1, _file2=None, _file3=None, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    file1 = os.path.join(proj_root, 'data/ent_N_flow', _model, _file1)

    with open(file1, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        ent_1 = []
        for i, row in enumerate(plots):
            ent_1.append(float(row[0]))

    if _file2 is None and _file3 is None:

        plt.figure(figsize=(5, 5))
        ax0 = plt.subplot(111)
        ax0.plot(ent_1, '.-', lw=1, label="$\delta=0.1$")
        ax0.set_xlabel('$N$')
        ax0.set_xlim([0, len(ent_1)])
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel('$S$')
        ax0.set_ylim([0, 0.6])
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('ent_N_flow_', '').replace('_', '\_').replace('.dat', ''))

        if _file3 is None and _save:
            os.makedirs(os.path.join(proj_root, 'figures/ent_N_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/ent_N_flow', _model, _file1.replace('.dat', '')+'.png'),
                        bbox_inches='tight', dpi=300)
        plt.show()

    elif _file3 is None:

        file2 = os.path.join(proj_root, 'data/ent_N_flow', _model, _file2)

        with open(file2, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            ent_2 = []
            for i, row in enumerate(plots):
                ent_2.append(float(row[0]))

        plt.figure(figsize=(10, 5))
        gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0.1)
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1], sharey=ax0)

        ax0.plot(ent_1, '.-', lw=1)
        ax0.set_xlabel('$N$')
        ax0.set_xlim([0, len(ent_1)])
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel('$S$')
        ax0.set_ylim([0, 0.6])
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('ent_N_flow_', '').replace('_', '\_').replace('.dat', ''))

        ax1.yaxis.set_visible(False)
        ax1.plot(ent_2, '.-', lw=1)
        ax1.set_xlabel('$N$')
        ax1.set_xlim([0, len(ent_2)])
        ax1.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax1.set_title(_file2.replace('ent_N_flow_', '').replace('_', '\_').replace('.dat', ''))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/ent_N_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/ent_N_flow', _model,
                                     _file1.replace('.dat', '')+'_comparison.png'), bbox_inches='tight', dpi=300)
        plt.show()

    else:

        file2 = os.path.join(proj_root, 'data/ent_N_flow', _model, _file2)

        with open(file2, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            ent_2 = []
            for i, row in enumerate(plots):
                ent_2.append(float(row[0]))

        file3 = os.path.join(proj_root, 'data/ent_N_flow', _model, _file3)

        with open(file3, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            ent_3 = []
            for i, row in enumerate(plots):
                ent_3.append(float(row[0]))

        plt.figure(figsize=(5, 5))
        ax0 = plt.subplot(111)
        ax0.plot(ent_1, '.-', c='C1', lw=1, label="$\delta=0.1$")
        ax0.plot(ent_2, '.-', c='C2', lw=1, label="$\delta=0.5$")
        ax0.plot(ent_3, '.-', c='C3', lw=1, label="$\delta=0.9$")
        ax0.set_xlabel('$N$')
        ax0.set_xlim([0, len(ent_1)])
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel('$S$')
        ax0.set_ylim([0, 0.6])
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('ent_N_flow_', '').replace('_', '\_').replace('.dat', ''))

        ax0.legend(loc='lower right')

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/ent_N_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/ent_N_flow', _model,
                                     _file1.replace('.dat', '')+'_comparison.png'), bbox_inches='tight', dpi=300)
        plt.show()


if __name__ == "__main__":

    model = 'spin2021'
    file1 = 'ent_N_flow_spin2021_L_8_obc_dis_10_J_1_1_1_T0_1_T1_1_N_1001_delta_0.9_W_2.dat'
    file2 = 'ent_N_flow_spin2021_L_8_obc_dis_10_J_1_1_1_T0_1_T1_1_N_1001_delta_0.1_W_2.dat'

    plot_ent_N_flow(model, file1, file2, _save=True)
