# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_ent_rel_W_flow(_model, _file1, _file2=None, _file3=None, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    with open(os.path.join(proj_root, 'data/ent_rel_W_flow', _model, _file1), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        W, ent_rel = [], []
        for i, row in enumerate(plots):
            W.append(float(row[0]))
            ent_rel.append(float(row[1]))

    plt.figure(figsize=(10, 5))
    ax0 = plt.subplot(111)
    ax0.plot(W, ent_rel, '.-', label="$L=8$")
    ax0.set_xlabel("$W$")
    # ax0.set_xscale('log', basex=2)
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel("$S_\mathrm{rel}$")
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_title(_file1.replace('ent_rel_W_flow_', '').replace('_', '\_').replace('.dat', ''))
    # ax0.axvline(x=3, ls='--', c='k')

    if _file2 is not None:

        with open(os.path.join(proj_root, 'data/ent_rel_W_flow', _model, _file2), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            W2, ent_rel2 = [], []
            for i, row in enumerate(plots):
                W2.append(float(row[0]))
                ent_rel2.append(float(row[1]))
        ax0.plot(W2, ent_rel2, '.-', label="$L=10$")

    if _file3 is not None:

        with open(os.path.join(proj_root, 'data/ent_rel_W_flow', _model, _file3), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            W3, ent_rel3 = [], []
            for i, row in enumerate(plots):
                W3.append(float(row[0]))
                ent_rel3.append(float(row[1]))
        ax0.plot(W3, ent_rel3, '.-', label="$L=12$")

    ax0.legend()

    if _save:
        os.makedirs(os.path.join(proj_root, 'figures/ent_rel_W_flow', _model), exist_ok=True)
        if _file2 is not None or _file3 is not None:
            plt.savefig(os.path.join(proj_root, 'figures/ent_rel_W_flow', _model,
                                     _file1.replace(".dat", "")+"_comparison.png"),
                        bbox_inches='tight', dpi=300)
        else:
            plt.savefig(os.path.join(proj_root, 'figures/ent_rel_W_flow', _model,
                                     _file1.replace(".dat", "")+".png"),
                        bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":

    model = 'heisenberg'
    file1 = 'ent_rel_W_flow_heisenberg_L_8_Nup_4_pauli_0_obc_dis_999_J_1_1_1_W_0.2_4.9_24.dat'
    file2 = 'ent_rel_W_flow_heisenberg_L_10_Nup_5_pauli_0_obc_dis_999_J_1_1_1_W_0.2_4.9_24.dat'
    # file3 = 'ent_rel_W_flow_heisenberg_L_12_Nup_6_pauli_0_obc_dis_1000_J_1_1_1_W_0.2_4.9_23.dat'

    plot_ent_rel_W_flow(model, file1, file2, _save=False)
