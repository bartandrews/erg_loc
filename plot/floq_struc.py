# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_floq_struc(_model, _file1, _file2=None, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    file1 = os.path.join(proj_root, 'data/floq_struc', _model, _file1)

    with open(file1, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        A2_1 = []
        for i, row in enumerate(plots):
            A2_1.append(float(row[0]))

    if _file2 is None:

        plt.figure(figsize=(5, 5))
        ax0 = plt.subplot(111)
        ax0.plot(A2_1, '.', c='k', lw=1)
        ax0.set_xlabel('$\\alpha$')
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel("$|A_{\\alpha, 0}|^2$")
        ax0.set_yscale('log')
        ax0.set_ylim([10e-15, 1])
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('floq_struc_', '').replace('_', '\_').replace('.dat', ''))
        ax0.axhline(0.0001, c='r')

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/floq_struc', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/floq_struc', _model, _file1.replace('.dat', '.png')),
                        bbox_inches='tight', dpi=300)
        plt.show()

    else:

        file2 = os.path.join(proj_root, 'data/floq_struc', _model, _file2)

        with open(file2, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            A2_2 = []
            for i, row in enumerate(plots):
                A2_2.append(float(row[0]))

        plt.figure(figsize=(10, 5))
        gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0)
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1], sharey=ax0)

        ax0.plot(A2_1, '.', c='k', lw=1)
        ax0.set_xlabel('$\\alpha$')
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel("$|A_{\\alpha, 0}|^2$")
        ax0.set_yscale('log')
        ax0.set_ylim([10e-15, 1])
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('floq_struc_', '').replace('_', '\_').replace('.dat', ''))
        ax0.axhline(0.0001, c='r')

        ax1.yaxis.set_visible(False)
        ax1.plot(A2_2, '.', c='k', lw=1)
        ax1.set_xlabel('$\\alpha$')
        ax1.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax1.set_title(_file2.replace('floq_struc_', '').replace('_', '\_').replace('.dat', ''))
        ax1.axhline(0.0001, c='r')

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/floq_struc', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/floq_struc', _model,
                                     _file2.replace('.dat', '_comparison.png')), bbox_inches='tight', dpi=300)
        plt.show()


if __name__ == "__main__":

    model = 'spin2021'
    file1 = 'floq_struc_spin2021_L_8_obc_J_1_1_1_T0_10_W_10.dat'
    # file2 = 'floq_struc_spin2021_L_8_obc_J_1_1_1_T0_10_W_8.dat'

    plot_floq_struc(model, file1, _save=False)
