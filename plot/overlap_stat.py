# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import os
from scipy.stats import gaussian_kde
import numpy as np
import scipy.stats as stats

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amsfonts}\usepackage{braket}')


def plot_overlap_stat(_model, _params1, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    L_list = [10, 12, 14, 16, 18, 20]

    plt.figure(figsize=(5, 5))
    ax0 = plt.subplot(111)

    for i, L_val in enumerate(L_list):

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

    ax0.legend(title="$L$")
    ax0.set_xlabel('$|\Delta Z_1|$')
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel('$p(|\Delta Z_1|)$')
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_title(f"{_model}_{_params1}".replace('_', '\_'))

    if _save:
        os.makedirs(os.path.join(proj_root, 'figures/overlap_stat', _model), exist_ok=True)
        plt.savefig(os.path.join(proj_root, 'figures/overlap_stat', _model, f"overlap_stat_{_model}_{_params1}.png"),
                    bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":

    model = 'pxp'
    params1 = 'obc_J_1_1_1_W_0'

    plot_overlap_stat(model, params1, _save=True)
