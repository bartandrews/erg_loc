# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import os
import numpy as np

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_ent_delta_flow(_model, _file1, _file2=None, _file3=None, _multi=False, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    with open(os.path.join(proj_root, 'data/ent_delta_flow', _model, _file1), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        nrow = len(list(plots))
        csvfile.seek(0)
        ncol = len(next(plots))-1
        csvfile.seek(0)
        delta = []
        ent = np.zeros((nrow, ncol), dtype=float)
        for i, row in enumerate(plots):
            delta.append(float(row[0]))
            for j in range(ncol):
                ent[i][j] = float(row[j+1])
        mean_ent = np.mean(ent, axis=1)

    plt.figure(figsize=(10, 5))
    ax0 = plt.subplot(111)
    if _multi:
        for j in range(ncol):
            ax0.plot(delta, ent[:, j], '-', lw=0.1)
        ax0.plot(delta, mean_ent, 'x-', c='k')
    else:
        ax0.plot(delta, mean_ent, '.-', label="$L=6$")
    ax0.set_xlabel("$\delta$")
    ax0.set_xlim([0, 1])
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel("$S$")
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_title(_file1.replace('delta_flow_', '').replace('_', '\_').replace('.dat', ''))

    if _file2 is not None:
        with open(os.path.join(proj_root, 'data/ent_delta_flow', _model, _file2), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            nrow = len(list(plots))
            csvfile.seek(0)
            ncol = len(next(plots)) - 1
            csvfile.seek(0)
            delta2 = []
            ent2 = np.zeros((nrow, ncol), dtype=float)
            for i, row in enumerate(plots):
                delta2.append(float(row[0]))
                for j in range(ncol):
                    ent2[i][j] = float(row[j + 1])
            mean_ent2 = np.mean(ent2, axis=1)
        ax0.plot(delta2, mean_ent2, '.-', label="$L=8$")
    if _file3 is not None:
        with open(os.path.join(proj_root, 'data/ent_delta_flow', _model, _file3), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            nrow = len(list(plots))
            csvfile.seek(0)
            ncol = len(next(plots)) - 1
            csvfile.seek(0)
            delta3 = []
            ent3 = np.zeros((nrow, ncol), dtype=float)
            for i, row in enumerate(plots):
                delta3.append(float(row[0]))
                for j in range(ncol):
                    ent3[i][j] = float(row[j + 1])
            mean_ent3 = np.mean(ent3, axis=1)
        ax0.plot(delta3, mean_ent3, '.-', label="$L=10$")

    if not _multi:
        ax0.legend(loc='upper left')

    if _save:
        if _file2 is None and _file3 is None:
            os.makedirs(os.path.join(proj_root, 'figures/ent_delta_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/ent_delta_flow', _model,
                                     _file1.replace(".dat", "")+".png"),
                        bbox_inches='tight', dpi=300)
        else:
            os.makedirs(os.path.join(proj_root, 'figures/ent_delta_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/ent_delta_flow', _model,
                                     _file1.replace(".dat", "") + "_comparison.png"),
                        bbox_inches='tight', dpi=300)

    plt.show()


if __name__ == "__main__":

    model = 'spin2021'
    file1 = 'ent_delta_flow_spin2021_L_6_obc_dis_100_J_1_1_1_T0_1_T1_1_delta_0_1_11_W_2.dat'
    file2 = 'ent_delta_flow_spin2021_L_8_obc_dis_100_J_1_1_1_T0_1_T1_1_delta_0_1_11_W_2.dat'
    file3 = 'ent_delta_flow_spin2021_L_10_obc_dis_100_J_1_1_1_T0_1_T1_1_delta_0_1_11_W_2.dat'

    plot_ent_delta_flow(model, file1, file2, file3, _multi=False, _save=True)
