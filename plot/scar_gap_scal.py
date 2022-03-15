# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import os
import numpy as np

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amsfonts}\usepackage{braket}')


def plot_scar_gap_scal(_model, _params1, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    L_list = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28]
    scar_ener = np.zeros((len(L_list), np.max(L_list)//2+1))
    scar_gaps = np.zeros((len(L_list), np.max(L_list)//2))
    gaps = np.zeros((len(L_list), 3))

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

        old_frac = 0
        for k, frac in enumerate(np.linspace(0.0415, 1.0415, L_val//2+1)):
            tmp_E_1_list, tmp_overlap_list = [], []
            for i, E_1_val in enumerate(E_1):
                if min(E_1)+old_frac*(max(E_1)-min(E_1)) <= E_1_val < min(E_1)+frac*(max(E_1)-min(E_1)):
                    tmp_E_1_list.append(E_1_val)
                    tmp_overlap_list.append(overlap_1[i])
            scar_ener[j, k] = tmp_E_1_list[np.argmax(tmp_overlap_list)]
            old_frac = frac

        for m in range(np.max(L_list)//2):
            if scar_ener[j, m+1] != 0:
                scar_gaps[j, m] = scar_ener[j, m+1] - scar_ener[j, m]

        for n in range(3):
            gap = scar_gaps[j, len([i for i in scar_gaps[j] if i != 0])//2 + n]
            gaps[j, n] = gap

    plt.figure(figsize=(10, 5))
    ax0 = plt.subplot(111)

    ax0.plot([1/i for i in L_list], gaps[:, 0], '.-', label="$\Delta E_1$")
    ax0.plot([1/i for i in L_list], gaps[:, 1], '.-', label="$\Delta E_2$")
    ax0.plot([1/i for i in L_list], gaps[:, 2], '.-', label="$\Delta E_3$")
    ax0.legend()
    ax0.set_xlabel('$1/L$')
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel('$\Delta E$')
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_title(f"{_model}_{_params1}".replace('_', '\_'))

    if _save:
        os.makedirs(os.path.join(proj_root, 'figures/scar_gap_scal', _model), exist_ok=True)
        plt.savefig(os.path.join(proj_root, 'figures/scar_gap_scal', _model, f"scar_gap_scal_{_model}_{_params1}.png"),
                    bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":

    model = 'pxp'
    params1 = 'obc_J_1_1_1_W_0'

    plot_scar_gap_scal(model, params1, _save=True)
