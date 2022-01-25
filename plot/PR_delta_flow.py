# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_PR_delta_flow(_model, _file1, _file2=None, _file3=None, _file4=None, _file5=None, _file6=None, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    with open(os.path.join(proj_root, 'data/PR_delta_flow', _model, _file1), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        delta, PR = [], []
        for i, row in enumerate(plots):
            delta.append(float(row[0]))
            PR.append(float(row[1]))

    plt.figure(figsize=(10, 5))
    ax0 = plt.subplot(111)
    ax0.plot(delta, PR, '.-', label="$L=100$")
    ax0.set_xlabel("$\delta$")
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel("PR")
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_title(_file1.replace('delta_flow_', '').replace('_', '\_').replace('.dat', ''))

    if _file2 is not None:
        with open(os.path.join(proj_root, 'data/PR_delta_flow', _model, _file2), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            delta2, PR2 = [], []
            for i, row in enumerate(plots):
                delta2.append(float(row[0]))
                PR2.append(float(row[1]))
        ax0.plot(delta2, PR2, '.-', label="$L=200$")
    if _file3 is not None:
        with open(os.path.join(proj_root, 'data/PR_delta_flow', _model, _file3), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            delta3, PR3 = [], []
            for i, row in enumerate(plots):
                delta3.append(float(row[0]))
                PR3.append(float(row[1]))
        ax0.plot(delta3, PR3, '.-', label="$L=300$")
    if _file4 is not None:
        with open(os.path.join(proj_root, 'data/PR_delta_flow', _model, _file4), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            delta4, PR4 = [], []
            for i, row in enumerate(plots):
                delta4.append(float(row[0]))
                PR4.append(float(row[1]))
        ax0.plot(delta4, PR4, '.-', label="$L=400$")
    if _file5 is not None:
        with open(os.path.join(proj_root, 'data/PR_delta_flow', _model, _file5), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            delta5, PR5 = [], []
            for i, row in enumerate(plots):
                delta5.append(float(row[0]))
                PR5.append(float(row[1]))
        ax0.plot(delta5, PR5, '.-', label="$L=500$")
    if _file6 is not None:
        with open(os.path.join(proj_root, 'data/PR_delta_flow', _model, _file6), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            delta6, PR6 = [], []
            for i, row in enumerate(plots):
                delta6.append(float(row[0]))
                PR6.append(float(row[1]))
        ax0.plot(delta6, PR6, '.-', label="$L=600$")

    ax0.legend(loc='upper left')

    if _save:
        if _file2 is None and _file3 is None:
            os.makedirs(os.path.join(proj_root, 'figures/PR_delta_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/PR_delta_flow', _model,
                                     _file1.replace(".dat", "")+".png"),
                        bbox_inches='tight', dpi=300)
        else:
            os.makedirs(os.path.join(proj_root, 'figures/PR_delta_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/PR_delta_flow', _model,
                                     _file1.replace(".dat", "") + "_comparison.png"),
                        bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":

    model = 'spin2021'
    file1 = 'PR_delta_flow_spin2021_L_100_Nup_1_obc_dis_10_J_1_1_1_T0_1_T1_1_delta_0_1_21_W_2.dat'
    file2 = 'PR_delta_flow_spin2021_L_200_Nup_1_obc_dis_10_J_1_1_1_T0_1_T1_1_delta_0_1_21_W_2.dat'
    file3 = 'PR_delta_flow_spin2021_L_300_Nup_1_obc_dis_10_J_1_1_1_T0_1_T1_1_delta_0_1_21_W_2.dat'
    file4 = 'PR_delta_flow_spin2021_L_400_Nup_1_obc_dis_10_J_1_1_1_T0_1_T1_1_delta_0_1_21_W_2.dat'
    file5 = 'PR_delta_flow_spin2021_L_500_Nup_1_obc_dis_10_J_1_1_1_T0_1_T1_1_delta_0_1_21_W_2.dat'
    file6 = 'PR_delta_flow_spin2021_L_600_Nup_1_obc_dis_10_J_1_1_1_T0_1_T1_1_delta_0_1_21_W_2.dat'

    plot_PR_delta_flow(model, file1, file2, file3, file4, file5, file6, _save=True)
