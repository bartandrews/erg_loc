# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def entropy_calculator(_proj_root, _model, _file):

    with open(os.path.join(_proj_root, 'data/ent_L_flow', _model, _file), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        L, S = [], []
        for i, row in enumerate(plots):
            L.append(float(row[0]))
            S.append(float(row[1]))
    dim = []
    for i, L_val in enumerate(L):
        if i == len(L)-1 or L_val != L[i+1]:
            dim.append(i+1)
    dim_cum = dim
    for i, dim_val in enumerate(dim):
        if i > 0:
            for j in range(i):
                dim[i] -= dim[j]
    S_ground, S_mid, S_max = [], [], []
    L_list = [int(i) for i in list(set(L))]
    for i in range(len(L_list)):
        if i == 0:
            S_ground.append(S[0])
            S_mid.append(S[int((-1+dim[0]/2))])
            S_max.append(S[dim_cum[i]-1])
        else:
            S_ground.append(S[dim_cum[i-1]])
            S_mid.append(S[int(dim_cum[i-1]-1+dim[i]/2)])
            S_max.append(S[dim_cum[i]-1])

    return L_list, S_ground, S_mid, S_max


def plot_ener_stat_W_flow(_model, _file1, _file2=None, _file3=None, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    plt.figure(figsize=(10, 5))
    ax0 = plt.subplot(111)
    (L_1_list, S_1_ground, S_1_mid, S_1_max) = entropy_calculator(proj_root, _model, _file1)
    ax0.plot(L_1_list, [x/y for x, y in zip(S_1_mid, L_1_list)], '.-', label="$W=0$ (mid)")
    if _file2 is not None:
        (L_2_list, S_2_ground, S_2_mid, S_2_max) = entropy_calculator(proj_root, _model, _file2)
        ax0.plot(L_2_list, [x / y for x, y in zip(S_2_mid, L_2_list)], '.-', label="$W=2.5$ (mid)")
    if _file3 is not None:
        (L_3_list, S_3_ground, S_3_mid, S_3_max) = entropy_calculator(proj_root, _model, _file3)
        ax0.plot(L_3_list, [x / y for x, y in zip(S_3_mid, L_3_list)], '.-', label="$W=5$ (mid)")
    ax0.legend()
    ax0.set_xlabel("$L$")
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel("$S/L$")
    ax0.set_yscale('log')
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_title(_file1.replace('ent_L_flow_', '').replace('_', '\_').replace('.dat', ''))

    if _save:
        os.makedirs(os.path.join(proj_root, 'figures/ent_L_flow', _model), exist_ok=True)
        if _file2 is not None or _file3 is not None:
            plt.savefig(os.path.join(proj_root, 'figures/ent_L_flow', _model,
                                     _file1.replace(".dat", "_comparison.png")),
                        bbox_inches='tight', dpi=300)
        else:
            plt.savefig(os.path.join(proj_root, 'figures/ent_L_flow', _model,
                                     _file1.replace(".dat", ".png")),
                        bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":

    model = 'heisenberg'
    file1 = 'ent_L_flow_heisenberg_L_8_12_3_obc_dis_10_J_1_1_1_W_0.dat'
    file2 = 'ent_L_flow_heisenberg_L_8_12_3_obc_dis_10_J_1_1_1_W_2.5.dat'
    file3 = 'ent_L_flow_heisenberg_L_8_12_3_obc_dis_10_J_1_1_1_W_5.dat'

    plot_ener_stat_W_flow(model, file1, file2, file3, _save=True)
