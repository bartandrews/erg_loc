# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_ent_L_flow(_model, _file1, _file2=None, _file3=None, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    plt.figure(figsize=(10, 5))
    ax0 = plt.subplot(111)
    with open(os.path.join(proj_root, 'data/ent_mid_L_flow', _model, _file1), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        L, S = [], []
        for i, row in enumerate(plots):
            L.append(float(row[0]))
            S.append(float(row[1]))
    ax0.plot(L, [y/x for x, y in zip(L, S)], '.-', label="$W=0$", c='tab:red')
    if _file2 is not None:
        with open(os.path.join(proj_root, 'data/ent_mid_L_flow', _model, _file2), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            L, S = [], []
            for i, row in enumerate(plots):
                L.append(float(row[0]))
                S.append(float(row[1]))
        ax0.plot(L, [y/x for x, y in zip(L, S)], '.-', label="$W=2.5$", c='tab:green')
    if _file3 is not None:
        with open(os.path.join(proj_root, 'data/ent_mid_L_flow', _model, _file3), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            L, S = [], []
            for i, row in enumerate(plots):
                L.append(float(row[0]))
                S.append(float(row[1]))
        ax0.plot(L, [y/x for x, y in zip(L, S)], '.-', label="$W=5$", c='tab:purple')
    ax0.legend()
    ax0.set_xlabel("$L$")
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel("$S/L$")
    ax0.set_yscale('log')
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_title(_file1.replace('ent_mid_L_flow_', '').replace('_', '\_').replace('.dat', ''))

    if _save:
        os.makedirs(os.path.join(proj_root, 'figures/ent_mid_L_flow', _model), exist_ok=True)
        if _file2 is not None or _file3 is not None:
            plt.savefig(os.path.join(proj_root, 'figures/ent_mid_L_flow', _model,
                                     _file1.replace(".dat", "_comparison.png")),
                        bbox_inches='tight', dpi=300)
        else:
            plt.savefig(os.path.join(proj_root, 'figures/ent_mid_L_flow', _model,
                                     _file1.replace(".dat", ".png")),
                        bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":

    model = 'heisenberg'
    file1 = 'ent_mid_L_flow_heisenberg_L_8_16_5_Nup_4_8_5_pauli_0_obc_dis_100_J_1_1_1_W_0.dat'
    file2 = 'ent_mid_L_flow_heisenberg_L_8_16_5_Nup_4_8_5_pauli_0_obc_dis_100_J_1_1_1_W_2.5.dat'
    file3 = 'ent_mid_L_flow_heisenberg_L_8_16_5_Nup_4_8_5_pauli_0_obc_dis_100_J_1_1_1_W_5.dat'

    plot_ent_L_flow(model, file1, file2, file3, _save=True)
