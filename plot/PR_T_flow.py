# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_PR_T_flow(_model, _file1, _file2=None, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    file1 = os.path.join(proj_root, 'data/PR_T_flow', _model, _file1)

    with open(file1, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        T_1, PR_1 = [], []
        for i, row in enumerate(plots):
            T_1.append(float(row[0]))
            PR_1.append(float(row[1]))

    if _file2 is None:

        plt.figure(figsize=(5, 5))
        ax0 = plt.subplot(111)
        ax0.plot(T_1, PR_1, '-', marker='x', c='k', lw=1)
        ax0.set_xlabel('$T$')
        ax0.set_xlim([0, max(T_1)])
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel('PR($i=0$)')
        # ax0.set_ylim([0, 1])
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('PR_T_flow_', '').replace('_', '\_').replace('.dat', ''))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/PR_T_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/PR_T_flow', _model, file1.replace('.dat', '.png')),
                        bbox_inches='tight', dpi=300)
        plt.show()

    else:

        file2 = os.path.join(proj_root, 'data/PR_T_flow', _model, _file2)

        with open(file2, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            T_2, PR_2 = [], []
            for i, row in enumerate(plots):
                T_2.append(float(row[0]))
                PR_2.append(float(row[1]))

        plt.figure(figsize=(10, 5))
        gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0.1)
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1], sharey=ax0)

        ax0.plot(T_1, PR_1, '.-', c='k', lw=1)
        ax0.set_xlabel('$T$')
        ax0.set_xlim([0, max(T_1)])
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel('PR($i=0$)')
        # ax0.set_ylim([0, 1])
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('PR_T_flow_', '').replace('_', '\_').replace('.dat', ''))

        ax1.yaxis.set_visible(False)
        ax1.plot(T_2, PR_2, '.-', c='k', lw=1)
        ax1.set_xlabel('$T$')
        ax1.set_xlim([0, max(T_2)])
        ax1.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax1.set_title(_file2.replace('PR_T_flow_', '').replace('_', '\_').replace('.dat', ''))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/PR_T_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/PR_T_flow', _model,
                                     _file1.replace('.dat', '_comparison.png')), bbox_inches='tight', dpi=300)
        plt.show()


if __name__ == "__main__":

    model = 'ponte2015'
    file1 = 'PR_T_flow_ponte2015_L_8_obc_dis_100_J_1_1_1_h0_2_T0_7_T_0_3_11_W_0.5.dat'
    file2 = 'PR_T_flow_ponte2015_L_8_obc_dis_100_J_1_1_1_h0_2_T0_7_T_0_3_11_W_8.dat'

    plot_PR_T_flow(model, file1, file2, _save=True)
