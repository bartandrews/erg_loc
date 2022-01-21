# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_r_W_flow(_model, _file, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    with open(os.path.join(proj_root, 'data/loc_len_delta_flow', _model, _file), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        delta, loc_len = [], []
        for i, row in enumerate(plots):
            delta.append(float(row[0]))
            loc_len.append(float(row[1]))

    plt.figure(figsize=(10, 5))
    ax0 = plt.subplot(111)
    ax0.plot(delta, loc_len, '.-')
    ax0.set_xlabel("$\delta$")
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel("$\\xi$")
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_title(_file.replace('delta_flow_', '').replace('_', '\_').replace('.dat', ''))

    if _save:
        os.makedirs(os.path.join(proj_root, 'figures/loc_len_delta_flow', _model), exist_ok=True)
        plt.savefig(os.path.join(proj_root, 'figures/loc_len_delta_flow', _model,
                                 _file.replace(".dat", "")+".png"),
                    bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":

    model = 'spin2021'
    file = 'loc_len_delta_flow_spin2021_L_200_Nup_1_obc_dis_10_J_1_1_1_T0_1_T1_1_delta_0_1_11_W_3.dat'

    plot_r_W_flow(model, file, _save=False)
