# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_ener_abs_N_flow(_model, _file1, _file2=None, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    file1 = os.path.join(proj_root, 'data/ener_abs_N_flow', _model, _file1)

    with open(file1, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        Q_1 = []
        for i, row in enumerate(plots):
            Q_1.append(float(row[0]))

    if _file2 is None:

        plt.figure(figsize=(5, 5))
        ax0 = plt.subplot(111)
        ax0.plot(Q_1, '-', marker='x', c='k', lw=1)
        ax0.set_xlabel('$N$')
        ax0.set_xlim([0, len(Q_1) - 1])
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel('$Q_N$')
        ax0.set_ylim([0, 1])
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('ener_abs_N_flow_', '').replace('_', '\_').replace('.dat', ''))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/ener_abs_N_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/ener_abs_N_flow', _model, file1.replace('.dat', '')+'.png'),
                        bbox_inches='tight', dpi=300)
        plt.show()

    else:

        file2 = os.path.join(proj_root, 'data/ener_abs_N_flow', _model, _file2)

        with open(file2, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            Q_2 = []
            for i, row in enumerate(plots):
                Q_2.append(float(row[0]))

        plt.figure(figsize=(10, 5))
        gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0.1)
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1], sharey=ax0)

        ax0.plot(Q_1, '.-', c='k', lw=1)
        ax0.set_xlabel('$N$')
        ax0.set_xlim([0, len(Q_1)-1])
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel('$Q_N$')
        ax0.set_ylim([0, 1])
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('ener_abs_N_flow_', '').replace('_', '\_').replace('.dat', ''))

        ax1.yaxis.set_visible(False)
        ax1.plot(Q_2, '.-', c='k', lw=1)
        ax1.set_xlabel('$N$')
        ax1.set_xlim([0, len(Q_2)-1])
        ax1.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax1.set_title(_file2.replace('ener_abs_N_flow_', '').replace('_', '\_').replace('.dat', ''))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/ener_abs_N_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/ener_abs_N_flow', _model,
                                     _file1.replace('.dat', '')+'_comparison.png'), bbox_inches='tight', dpi=300)
        plt.show()


if __name__ == "__main__":

    model = 'ponte2015'
    file1 = 'ener_abs_N_flow_ponte2015_L_8_obc_dis_100_J_1_1_1_h0_2_T0_7_T1_1.5_N_31_W_0.5.dat.new'
    file2 = 'ener_abs_N_flow_ponte2015_L_8_obc_dis_100_J_1_1_1_h0_2_T0_7_T1_1.5_N_31_W_8.dat.new'

    plot_ener_abs_N_flow(model, file1, file2, _save=False)
