# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_q_ener_spec(_model, _file1, _file2=None, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    with open(os.path.join(proj_root, 'data/q_ener', _model, _file1), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        E_1 = []
        for i, row in enumerate(plots):
            E_1.append(float(row[0]))

    if _file2 is None:

        plt.figure(figsize=(5, 5))
        ax0 = plt.subplot(111)
        ax0.plot([0]*len(E_1), E_1, '.', marker='_', c='k', lw=1)
        ax0.xaxis.set_visible(False)
        ax0.set_ylabel('$E$')
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('q_ener_', '').replace('_', '\_').replace('.dat', ''))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/q_ener_spec', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/q_ener_spec', _model,
                                     _file1.replace("q_ener_", "q_ener_spec_").replace(".dat", "")+".png"),
                        bbox_inches='tight', dpi=300)
        plt.show()

    else:

        with open(os.path.join(proj_root, 'data/q_ener', _model, _file2), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            E_2 = []
            for i, row in enumerate(plots):
                E_2.append(float(row[0]))

        plt.figure(figsize=(10, 5))
        gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0)
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1], sharey=ax0)

        ax0.plot([0] * len(E_1), E_1, '.', marker='_', c='k', lw=1)
        ax0.xaxis.set_visible(False)
        ax0.set_ylabel('$E$')
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('q_ener_', '').replace('_', '\_').replace('.dat', ''))

        ax1.yaxis.set_visible(False)
        ax1.plot([0] * len(E_2), E_2, '.', marker='_', c='k', lw=1)
        ax1.xaxis.set_visible(False)
        ax1.set_ylabel('$E$')
        ax1.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax1.set_title(_file2.replace('q_ener_', '').replace('_', '\_').replace('.dat', ''))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/q_ener_spec', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/q_ener_spec', _model,
                                     _file1.replace("q_ener_", "q_ener_spec_").replace(".dat", "")+"_comparison.png"),
                        bbox_inches='tight', dpi=300)
        plt.show()


if __name__ == "__main__":

    model = 'spin2021'
    file1 = 'q_ener_spin2021_L_200_Nup_1_obc_dis_100_J_1_1_1_T0_1_T1_1_delta_1_W_3.dat'
    # file2 = 'q_ener_ponte2015_L_8_pauli_0_obc_J_1_1_1_h0_2_T0_7_T1_1.5_W_8.dat'

    plot_q_ener_spec(model, file1, _save=True)
