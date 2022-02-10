# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_r_W_flow(_model, _file, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    with open(os.path.join(proj_root, 'data/r_W_flow', _model, _file), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        W, r = [], []
        for i, row in enumerate(plots):
            W.append(float(row[0]))
            r.append(float(row[1]))

    plt.figure(figsize=(10, 5))
    ax0 = plt.subplot(111)
    ax0.plot(W, r, '.-')
    ax0.plot([min(W), max(W)], [0.53, 0.53], label="GOE", c='g')
    ax0.plot([min(W), max(W)], [0.39, 0.39], label="Poisson", c='r')
    ax0.legend()
    ax0.set_xlabel("$W$")
    # ax0.set_xscale('log', basex=2)
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel("$[r^{(n)}_\\alpha]$")
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_title(_file.replace('r_W_flow_', '').replace('_', '\_').replace('.dat', ''))
    # ax0.axvline(x=3, ls='--', c='k')

    if _save:
        os.makedirs(os.path.join(proj_root, 'figures/r_W_flow', _model), exist_ok=True)
        plt.savefig(os.path.join(proj_root, 'figures/r_W_flow', _model,
                                 _file.replace(".dat", "")+".png"),
                    bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":

    model = 'heisenberg'
    file = 'r_W_flow_heisenberg_L_8_Nup_4_pauli_0_obc_dis_1000_J_1_1_1_W_0.2_4.9_23.dat'

    plot_r_W_flow(model, file, _save=False)
