# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_ent_L_flow(_model, _file1, _file2=None, _file3=None, _file4=None, _file5=None, _file6=None,
                    _file7=None, _file8=None, _file9=None, _file10=None, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    plt.figure(figsize=(10, 5))
    ax0 = plt.subplot(111)
    with open(os.path.join(proj_root, 'data/ent_mid_L_flow', _model, _file1), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        L, S, S_page = [], [], []
        for i, row in enumerate(plots):
            L.append(float(row[0]))
            S.append(float(row[1]))
            S_page.append(float(row[2]))
    ax0.plot(L, [s for s, l in zip(S_page, L)], 'x-', label="Page", c='k')
    ax0.plot(L, [s for s, l in zip(S, L)], '.-', label="$W=0.5$", c='C0')
    if _file2 is not None:
        with open(os.path.join(proj_root, 'data/ent_mid_L_flow', _model, _file2), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            L, S = [], []
            for i, row in enumerate(plots):
                L.append(float(row[0]))
                S.append(float(row[1]))
        ax0.plot(L, [s for s, l in zip(S, L)], '.-', label="$W=1$", c='C1')
    if _file3 is not None:
        with open(os.path.join(proj_root, 'data/ent_mid_L_flow', _model, _file3), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            L, S = [], []
            for i, row in enumerate(plots):
                L.append(float(row[0]))
                S.append(float(row[1]))
        ax0.plot(L, [s for s, l in zip(S, L)], '.-', label="$W=1.5$", c='C2')
    if _file4 is not None:
        with open(os.path.join(proj_root, 'data/ent_mid_L_flow', _model, _file4), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            L, S = [], []
            for i, row in enumerate(plots):
                L.append(float(row[0]))
                S.append(float(row[1]))
        ax0.plot(L, [s for s, l in zip(S, L)], '.-', label="$W=2$", c='C3')
    if _file5 is not None:
        with open(os.path.join(proj_root, 'data/ent_mid_L_flow', _model, _file5), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            L, S = [], []
            for i, row in enumerate(plots):
                L.append(float(row[0]))
                S.append(float(row[1]))
        ax0.plot(L, [s for s, l in zip(S, L)], '.-', label="$W=2.5$", c='C4')
    if _file6 is not None:
        with open(os.path.join(proj_root, 'data/ent_mid_L_flow', _model, _file6), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            L, S = [], []
            for i, row in enumerate(plots):
                L.append(float(row[0]))
                S.append(float(row[1]))
        ax0.plot(L, [s for s, l in zip(S, L)], '.-', label="$W=3$", c='C5')
    if _file7 is not None:
        with open(os.path.join(proj_root, 'data/ent_mid_L_flow', _model, _file7), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            L, S = [], []
            for i, row in enumerate(plots):
                L.append(float(row[0]))
                S.append(float(row[1]))
        ax0.plot(L, [s for s, l in zip(S, L)], '.-', label="$W=3.5$", c='C6')
    if _file8 is not None:
        with open(os.path.join(proj_root, 'data/ent_mid_L_flow', _model, _file8), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            L, S = [], []
            for i, row in enumerate(plots):
                L.append(float(row[0]))
                S.append(float(row[1]))
        ax0.plot(L, [s for s, l in zip(S, L)], '.-', label="$W=4$", c='C7')
    if _file9 is not None:
        with open(os.path.join(proj_root, 'data/ent_mid_L_flow', _model, _file9), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            L, S = [], []
            for i, row in enumerate(plots):
                L.append(float(row[0]))
                S.append(float(row[1]))
        ax0.plot(L, [s for s, l in zip(S, L)], '.-', label="$W=4.5$", c='C8')
    if _file10 is not None:
        with open(os.path.join(proj_root, 'data/ent_mid_L_flow', _model, _file10), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            L, S = [], []
            for i, row in enumerate(plots):
                L.append(float(row[0]))
                S.append(float(row[1]))
        ax0.plot(L, [s for s, l in zip(S, L)], '.-', label="$W=5$", c='C9')
    ax0.legend(loc='best')
    ax0.set_xlabel("$L$")
    # ax0.set_xscale('log', basex=2)
    ax0.set_xlim([L[0], L[-1]])
    ax0.tick_params('x', direction='in', bottom=True, top=True)
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel("$S$")
    # ax0.set_yscale('log')
    ax0.tick_params('y', direction='in', left=True, right=True)
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_title(_file1.replace('ent_mid_L_flow_', '').replace('_', '\_').replace('.dat', ''))

    if _save:
        os.makedirs(os.path.join(proj_root, 'figures/ent_mid_L_flow', _model), exist_ok=True)
        if _file2 is not None or _file3 is not None or _file4 is not None or _file5 is not None or _file6 is not None\
                or _file7 is not None or _file8 is not None or _file9 is not None or _file10 is not None:
            plt.savefig(os.path.join(proj_root, 'figures/ent_mid_L_flow', _model,
                                     _file1.replace(".dat", "")+"_comparison.png"),
                        bbox_inches='tight', dpi=300)
        else:
            plt.savefig(os.path.join(proj_root, 'figures/ent_mid_L_flow', _model,
                                     _file1.replace(".dat", "")+".png"),
                        bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":

    model = 'heisenberg'
    file1 = 'ent_mid_L_flow_heisenberg_L_12_16_3_Nup_6_8_3_pauli_0_obc_dis_100_J_1_1_1_W_0.5.dat'
    file2 = 'ent_mid_L_flow_heisenberg_L_12_16_3_Nup_6_8_3_pauli_0_obc_dis_100_J_1_1_1_W_1.dat'
    file3 = 'ent_mid_L_flow_heisenberg_L_12_16_3_Nup_6_8_3_pauli_0_obc_dis_100_J_1_1_1_W_1.5.dat'
    file4 = 'ent_mid_L_flow_heisenberg_L_12_16_3_Nup_6_8_3_pauli_0_obc_dis_100_J_1_1_1_W_2.dat'
    file5 = 'ent_mid_L_flow_heisenberg_L_12_16_3_Nup_6_8_3_pauli_0_obc_dis_100_J_1_1_1_W_2.5.dat'
    file6 = 'ent_mid_L_flow_heisenberg_L_12_16_3_Nup_6_8_3_pauli_0_obc_dis_100_J_1_1_1_W_3.dat'
    file7 = 'ent_mid_L_flow_heisenberg_L_12_16_3_Nup_6_8_3_pauli_0_obc_dis_100_J_1_1_1_W_3.5.dat'
    file8 = 'ent_mid_L_flow_heisenberg_L_12_16_3_Nup_6_8_3_pauli_0_obc_dis_100_J_1_1_1_W_4.dat'
    file9 = 'ent_mid_L_flow_heisenberg_L_12_16_3_Nup_6_8_3_pauli_0_obc_dis_100_J_1_1_1_W_4.5.dat'
    file10 = 'ent_mid_L_flow_heisenberg_L_12_16_3_Nup_6_8_3_pauli_0_obc_dis_100_J_1_1_1_W_5.dat'

    plot_ent_L_flow(model, file1, file2, file3, file4, file5, file6, file7, file8, file9, file10, _save=True)
    # plot_ent_L_flow(model, file1, _save=True)
