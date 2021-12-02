# --- python imports
import matplotlib.pyplot as plt
import csv
import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_q_ener_stat(_model, _file1, _file2=None, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    with open(os.path.join(proj_root, 'data/q_ener_spac', _model, _file1), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        E_spac_1 = []
        for i, row in enumerate(plots):
            E_spac_1.append(float(row[0]))

    r_1 = []
    for i in range(1, len(E_spac_1)):
        r_1.append(E_spac_1[i]/E_spac_1[i-1])

    if _file2 is None:

        plt.figure(figsize=(5, 5))
        ax0 = plt.subplot(111)
        ax0.hist(r_1, bins=np.arange(0, 5.1, 0.1), density=True)
        r_vals = np.arange(0, 100, 0.1)
        Poisson = [1 / ((1 + r) ** 2) for r in r_vals]
        GOE = [(27 / 8) * ((r + r ** 2) / (1 + r + r ** 2) ** (5 / 2)) for r in r_vals]
        ax0.plot(r_vals, Poisson, c='r', label='Poisson')
        ax0.plot(r_vals, GOE, c='g', label='GOE')
        ax0.legend(loc='upper right')
        ax0.set_xlabel("$r$")
        ax0.set_xlim([0, 5])
        ax0.set_ylabel("$P(r)$")
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('q_ener_spac_', '').replace('_', '\_').replace('.dat', ''))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/q_ener_stat', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/q_ener_stat', _model,
                                     _file1.replace('spac', 'stat').replace(".dat", ".png")),
                        bbox_inches='tight', dpi=300)
        plt.show()

    else:

        with open(os.path.join(proj_root, 'data/q_ener_spac', _model, _file2), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            E_spac_2 = []
            for i, row in enumerate(plots):
                E_spac_2.append(float(row[0]))

        r_2 = []
        for i in range(1, len(E_spac_2)):
            r_2.append(E_spac_2[i]/E_spac_2[i-1])

        plt.figure(figsize=(10, 5))
        gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0.1)
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1], sharey=ax0)

        ax0.hist(r_1, bins=np.arange(0, 5.1, 0.1), density=True)
        r_vals = np.arange(0, 100, 0.1)
        Poisson = [1 / ((1 + r) ** 2) for r in r_vals]
        GOE = [(27 / 8) * ((r + r ** 2) / (1 + r + r ** 2) ** (5 / 2)) for r in r_vals]
        ax0.plot(r_vals, Poisson, c='r', label='Poisson')
        ax0.plot(r_vals, GOE, c='g', label='GOE')
        ax0.legend(loc='upper right')
        ax0.set_xlabel("$r$")
        ax0.set_xlim([0, 5])
        ax0.set_ylabel("$P(r)$")
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('q_ener_spac_', '').replace('_', '\_').replace('.dat', ''))

        ax1.yaxis.set_visible(False)
        ax1.hist(r_2, bins=np.arange(0, 5.1, 0.1), density=True)
        ax1.plot(r_vals, Poisson, c='r', label='Poisson')
        ax1.plot(r_vals, GOE, c='g', label='GOE')
        ax1.legend(loc='upper right')
        ax1.set_xlabel("$r$")
        ax1.set_xlim([0, 5])
        ax1.set_title(_file2.replace('q_ener_spac_', '').replace('_', '\_').replace('.dat', ''))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/q_ener_stat', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/q_ener_stat', _model,
                                     _file1.replace('spac', 'stat').replace(".dat", "_comparison.png")),
                        bbox_inches='tight', dpi=300)
        plt.show()


if __name__ == "__main__":

    model = 'spin2021'
    file1 = 'q_ener_spac_spin2021_L_8_Nup_4_obc_dis_100_J_1_1_1_T0_10_W_0.5.dat'
    file2 = 'q_ener_spac_spin2021_L_8_Nup_4_obc_dis_100_J_1_1_1_T0_10_W_8.dat'

    plot_q_ener_stat(model, file1, file2, _save=True)
