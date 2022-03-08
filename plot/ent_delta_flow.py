# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import os
import numpy as np
import h5py
import glob
import re
# --- QuSpin imports
from quspin.basis import spin_basis_1d

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_ent_delta_flow(_model, _file1, _file2=None, _file3=None, _file4=None, _multi=False, _save=False,
                        _h5_file=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    if _h5_file:

        L_pattern = re.search('L_(\d+)', _file1)
        L = int(L_pattern.group(1))

        basis = spin_basis_1d(L)

        delta_pattern = re.search('delta_(\d+)_(\d+)_(\d+)', _file1)
        delta_string = delta_pattern.group(0)
        delta_min = float(delta_pattern.group(1))
        delta_max = float(delta_pattern.group(2))
        delta_samp = int(delta_pattern.group(3))

        dis_pattern = re.search('dis_(\d+)', _file1)
        dis = int(dis_pattern.group(1))

        delta = []
        mean_ent = []
        for delta_val in np.linspace(delta_min, delta_max, delta_samp):  # delta

            delta.append(delta_val)

            _file1_tmp = _file1.replace("ent_delta_flow", "eig_U").replace(".dat", ".h5")\
                .replace(delta_string, f"delta_{delta_val:g}")
            _path1 = os.path.join(proj_root, 'hdf5/eig_U', _model, _file1_tmp)
            _path1_list = glob.glob(_path1.replace("_J", "_bat_*_J"))

            # determine number of states
            hf = h5py.File(_path1_list[0], 'r')
            U_array = hf['eig_U']
            Ns = np.shape(U_array)[2]

            _ent_array = np.zeros((dis, len(_path1_list), Ns))

            for dis_val in range(dis):  # disorder
                for i, _path1_list_val in enumerate(_path1_list):  # batch
                    hf = h5py.File(_path1_list_val, 'r')
                    U_array = hf['eig_U']
                    U_array = np.delete(U_array, 0, axis=1)
                    for k in range(Ns):  # state
                        _ent_array[dis_val, i, k] = basis.ent_entropy(U_array[dis_val, :, k],
                                                                      sub_sys_A=range(basis.L//2),
                                                                      density=False)["Sent_A"]
            mean_ent.append(np.mean(_ent_array, axis=(0, 1, 2)))
            print(mean_ent)
    else:
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
        for j, delta_val in enumerate(delta):
            ax0.plot([delta_val]*np.shape(ent)[1], ent[j, :], '_', markeredgewidth=0.2, c=f"C{j}")
        ax0.plot(delta, mean_ent, 'x-', c='k')
    else:
        for j, delta_val in enumerate(delta):
            ax0.errorbar(delta_val, np.mean(ent[j, :]), yerr=np.std(ent[j, :]), capsize=2, c="C0")
        ax0.plot(delta, mean_ent, '.-', label="$L=6$", c="C0")
    ax0.set_xlabel("$\delta$")
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
        for j, delta2_val in enumerate(delta2):
            ax0.errorbar(delta2_val, np.mean(ent2[j, :]), yerr=np.std(ent2[j, :]), capsize=2, c="C1")
        ax0.plot(delta2, mean_ent2, '.-', label="$L=8$", c="C1")
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
        for j, delta3_val in enumerate(delta3):
            ax0.errorbar(delta3_val, np.mean(ent3[j, :]), yerr=np.std(ent3[j, :]), capsize=2, c="C2")
        ax0.plot(delta3, mean_ent3, '.-', label="$L=10$", c="C2")
    if _file4 is not None:
        with open(os.path.join(proj_root, 'data/ent_delta_flow', _model, _file4), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            nrow = len(list(plots))
            csvfile.seek(0)
            ncol = len(next(plots)) - 1
            csvfile.seek(0)
            delta4 = []
            ent4 = np.zeros((nrow, ncol), dtype=float)
            for i, row in enumerate(plots):
                delta4.append(float(row[0]))
                for j in range(ncol):
                    ent4[i][j] = float(row[j + 1])
            mean_ent4 = np.mean(ent4, axis=1)
        for j, delta4_val in enumerate(delta4):
            ax0.errorbar(delta4_val, np.mean(ent4[j, :]), yerr=np.std(ent4[j, :]), capsize=2, c="C3")
        ax0.plot(delta4, mean_ent4, '.-', label="$L=12$", c="C3")

    if not _multi:
        ax0.legend(loc='upper left')

    if _save:
        if _file2 is None and _file3 is None and _file4 is None:
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
    file4 = 'ent_delta_flow_spin2021_L_12_obc_dis_100_J_1_1_1_T0_1_T1_1_delta_0_1_11_W_2.dat'

    plot_ent_delta_flow(model, file1, file2, file3, file4, _multi=False, _save=True, _h5_file=False)
