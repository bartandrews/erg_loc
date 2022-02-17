# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import os
import numpy as np

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_ent_L_flow(_model, _file1, _file2=None, _file3=None, _file4=None, _multi=False, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    with open(os.path.join(proj_root, 'data/ent_L_flow', _model, _file1), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        nrow = len(list(plots))
        csvfile.seek(0)
        ncol = len(next(plots))-1
        csvfile.seek(0)
        L = []
        ent = np.zeros((nrow, ncol), dtype=float)
        for i, row in enumerate(plots):
            L.append(float(row[0]))
            for j in range(ncol):
                ent[i][j] = float(row[j+1])
        mean_ent = np.mean(ent, axis=1)

    plt.figure(figsize=(10, 5))
    ax0 = plt.subplot(111)
    if _multi:
        for j in range(ncol):
            ax0.plot(L, ent[:, j], '-', lw=0.1)
        ax0.plot(L, mean_ent, 'x-', c='k')
    else:
        ax0.plot(L, mean_ent, '.-', label="$\delta=0.1$")
    ax0.set_xlabel("$L$")
    ax0.set_xlim([min(L), max(L)])
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel("$S$")
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_title(_file1.replace('L_flow_', '').replace('_', '\_').replace('.dat', ''))

    if _file2 is not None:
        with open(os.path.join(proj_root, 'data/ent_L_flow', _model, _file2), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            nrow = len(list(plots))
            csvfile.seek(0)
            ncol = len(next(plots)) - 1
            csvfile.seek(0)
            L2 = []
            ent2 = np.zeros((nrow, ncol), dtype=float)
            for i, row in enumerate(plots):
                L2.append(float(row[0]))
                for j in range(ncol):
                    ent2[i][j] = float(row[j + 1])
            mean_ent2 = np.mean(ent2, axis=1)
        ax0.plot(L2, mean_ent2, '.-', label="$\delta=0.2$")
    if _file3 is not None:
        with open(os.path.join(proj_root, 'data/ent_L_flow', _model, _file3), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            nrow = len(list(plots))
            csvfile.seek(0)
            ncol = len(next(plots)) - 1
            csvfile.seek(0)
            L3 = []
            ent3 = np.zeros((nrow, ncol), dtype=float)
            for i, row in enumerate(plots):
                L3.append(float(row[0]))
                for j in range(ncol):
                    ent3[i][j] = float(row[j + 1])
            mean_ent3 = np.mean(ent3, axis=1)
        ax0.plot(L3, mean_ent3, '.-', label="$\delta=0.3$")
    if _file4 is not None:
        with open(os.path.join(proj_root, 'data/ent_L_flow', _model, _file4), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            nrow = len(list(plots))
            csvfile.seek(0)
            ncol = len(next(plots)) - 1
            csvfile.seek(0)
            L4 = []
            ent4 = np.zeros((nrow, ncol), dtype=float)
            for i, row in enumerate(plots):
                L4.append(float(row[0]))
                for j in range(ncol):
                    ent4[i][j] = float(row[j + 1])
            mean_ent4 = np.mean(ent4, axis=1)
        ax0.plot(L4, mean_ent4, '.-', label="$\delta=0.4$")

    if not _multi:
        ax0.legend(loc='upper left')

    if _save:
        if _file2 is None and _file3 is None and _file4 is not None:
            os.makedirs(os.path.join(proj_root, 'figures/ent_L_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/ent_L_flow', _model,
                                     _file1.replace(".dat", "")+".png"),
                        bbox_inches='tight', dpi=300)
        else:
            os.makedirs(os.path.join(proj_root, 'figures/ent_L_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/ent_L_flow', _model,
                                     _file1.replace(".dat", "") + "_comparison.png"),
                        bbox_inches='tight', dpi=300)

    plt.show()


if __name__ == "__main__":

    model = 'spin2021'
    file1 = 'ent_L_flow_spin2021_L_8_24_5_Nup_4_12_5_obc_dis_10_J_1_1_1_T0_1_T1_1_delta_0.1_W_2.dat'
    # file2 = 'ent_L_flow_spin2021_L_8_obc_dis_100_J_1_1_1_T0_1_T1_1_delta_0_1_11_W_2.dat'
    # file3 = 'ent_L_flow_spin2021_L_10_obc_dis_100_J_1_1_1_T0_1_T1_1_delta_0_1_11_W_2.dat'
    # file4 = 'ent_L_flow_spin2021_L_12_obc_dis_100_J_1_1_1_T0_1_T1_1_delta_0_1_11_W_2.dat'

    plot_ent_L_flow(model, file1, _multi=True, _save=False)
