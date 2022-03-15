# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import os
from scipy.stats import gaussian_kde
import numpy as np
import scipy.stats as stats
# --- QuSpin imports
import functions.func_ham as fh

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amsfonts}\usepackage{braket}')


def plot_overlap_stat(_model, _params1, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    L_list = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28]

    plt.figure(figsize=(10, 5))
    gs = gridspec.GridSpec(1, 2, hspace=0)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    insert = np.zeros((len(L_list), 2))

    for j, L_val in enumerate(L_list):

        ener_file1 = os.path.join(proj_root, 'data/ener', _model, f"ener_{_model}_L_{L_val}_{_params1}.dat")
        overlap_file1 = os.path.join(proj_root, 'data/overlap', _model, f"overlap_{_model}_L_{L_val}_{_params1}.dat")

        if os.path.exists(ener_file1) is False or os.path.exists(overlap_file1) is False:
            raise ValueError("ener_file1 or overlap_file1 do not exist")

        with open(ener_file1, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            E_1 = []
            for i, row in enumerate(plots):
                E_1.append(float(row[0]))
        with open(overlap_file1, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            overlap_1 = []
            for i, row in enumerate(plots):
                overlap_1.append(float(row[0]))

        overlap_spac_1 = []
        for i, overlap_val in enumerate(overlap_1[:-1]):
            overlap_spac_1.append(np.abs(overlap_1[i+1]-overlap_1[i]))

        density_1 = stats.gaussian_kde(overlap_spac_1)

        # n, x, _ = ax0.hist(overlap_spac_1, bins=np.arange(0, 0.1, 0.005), density=True, histtype='step')
        x = np.arange(0, 0.15, 0.005)
        ax0.plot(x, density_1(x), label=f"${L_val}$")

        leaf_args = {"L": L_val, "Nup": None, "pauli": 0, "bc": "o", "J": (0, 0, 0)}
        Ns = fh.chosen_hamiltonian("pxp", leaf_args).Ns
        insert[j, 0] = Ns
        insert[j, 1] = np.mean(overlap_spac_1)

    ax0.legend(title="$L$")
    ax0.set_xlabel('$|\Delta Z_1|$')
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel('$p(|\Delta Z_1|)$')
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_title(f"{_model}_{_params1}".replace('_', '\_'))

    x, y = np.log(insert[:, 0]), np.log(insert[:, 1])
    ax1.plot(x, y, '.-')
    theta = np.polyfit(np.log(insert[:, 0]), np.log(insert[:, 1]), 1)
    y_line = theta[1] + theta[0] * x
    ax1.plot(x, y_line, c='r')
    ax1.set_xlabel('$\ln d_L$')
    ax1.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax1.set_ylabel('$\ln \overline{|\Delta Z_1|}$')
    ax1.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax1.set_title(f"$\overline{{|\Delta Z_1|}} \sim d_L^{{{theta[0]:.3g}}}$")

    if _save:
        os.makedirs(os.path.join(proj_root, 'figures/overlap_stat', _model), exist_ok=True)
        plt.savefig(os.path.join(proj_root, 'figures/overlap_stat', _model, f"overlap_stat_{_model}_{_params1}.png"),
                    bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":

    model = 'pxp'
    params1 = 'obc_J_1_1_1_W_0'

    plot_overlap_stat(model, params1, _save=True)
