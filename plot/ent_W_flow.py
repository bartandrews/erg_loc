# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_ent_W_flow(_model, _file, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    with open(os.path.join(proj_root, 'data/ent_W_flow', _model, _file), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        W, S = [], []
        for i, row in enumerate(plots):
            W.append(float(row[0]))
            S.append(float(row[1]))

    dim = 0
    for i, W_val in enumerate(W):
        if W_val != W[i+1]:
            dim = i+1
            break

    plt.figure(figsize=(10, 5))
    ax0 = plt.subplot(111)
    ax0.plot(W[0::dim], S[0::dim], '.-', label="$S_\mathrm{min}$")
    ax0.plot(W[int((dim/2)-1)::dim], S[int((dim/2)-1)::dim], '.-', label="$S_\mathrm{mid}$")
    ax0.plot(W[dim-1::dim], S[dim-1::dim], '.-', label="$S_\mathrm{max}$")
    ax0.legend()
    ax0.set_xlabel("$W$")
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel("$S$")
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_title(_file.replace('ent_W_flow_', '').replace('_', '\_').replace('.dat', ''))
    ax0.axvline(x=3, ls='--', c='k')

    if _save:
        os.makedirs(os.path.join(proj_root, 'figures/ent_W_flow', _model), exist_ok=True)
        plt.savefig(os.path.join(proj_root, 'figures/ent_W_flow', _model,
                                 _file.replace(".dat", "")+".png"),
                    bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":

    model = 'heisenberg'
    file = 'ent_W_flow_heisenberg_L_8_Nup_4_pauli_0_obc_dis_10000_J_1_1_1_W_0.5_12.5_24.dat'

    plot_ent_W_flow(model, file, _save=True)
