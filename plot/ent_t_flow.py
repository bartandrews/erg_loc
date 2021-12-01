# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_ent_t_flow(_model, _file, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    with open(os.path.join(proj_root, 'data/ent_t_flow', _model, _file), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        t, S = [], []
        for i, row in enumerate(plots):
            if i != 0:  # skip t=0 point
                t.append(float(row[0]))
                S.append(float(row[1]))

    plt.figure(figsize=(10, 5))
    ax0 = plt.subplot(111)
    ax0.plot(t, S, '.-', label="mid state")
    ax0.set_xlabel("$t$")
    ax0.set_xscale('log')
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel("$S$")
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_title(_file.replace('ent_t_flow_', '').replace('_', '\_').replace('.dat', ''))

    if _save:
        os.makedirs(os.path.join(proj_root, 'figures/ent_t_flow', _model), exist_ok=True)
        plt.savefig(os.path.join(proj_root, 'figures/ent_t_flow', _model,
                                 _file.replace(".dat", ".png")),
                    bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":

    model = 'heisenberg'
    file = 'ent_t_flow_heisenberg_L_6_obc_dis_100_t_-1_3_100_J_1_1_0.2_W_10.dat'

    plot_ent_t_flow(model, file, _save=True)
