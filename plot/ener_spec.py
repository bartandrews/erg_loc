# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_ener_spec(_model, _file1, _file2=None, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    with open(os.path.join(proj_root, 'data/ener_spec', _model, _file1), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        E_1 = []
        for i, row in enumerate(plots):
            E_1.append(float(row[0]))

    if _file2 is None:

        plt.figure(figsize=(5, 5))
        ax0 = plt.subplot(111)
        ax0.plot([0]*len(E_1), E_1, '.', marker='_', c='k', lw=1)
        ax0.xaxis.set_visible(False)
        ax0.set_ylabel('$E$')
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('ener_spec_', '').replace('_', '\_').replace('.dat', ''))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/ener_spec', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/ener_spec', _model, _file1.replace(".dat", "")+".png"),
                        bbox_inches='tight', dpi=300)
        plt.show()

    else:

        with open(os.path.join(proj_root, 'data/ener_spec', _model, _file2), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            E_2 = []
            for i, row in enumerate(plots):
                E_2.append(float(row[0]))

        plt.figure(figsize=(10, 5))
        gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0)
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1], sharey=ax0)

        ax0.plot([0] * len(E_1), E_1, '.', marker='_', c='k', lw=1)
        ax0.xaxis.set_visible(False)
        ax0.set_ylabel('$E$')
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('ener_spec_', '').replace('_', '\_').replace('.dat', ''))

        ax1.yaxis.set_visible(False)
        ax1.plot([0] * len(E_2), E_2, '.', marker='_', c='k', lw=1)
        ax1.xaxis.set_visible(False)
        ax1.set_ylabel('$E$')
        ax1.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax1.set_title(_file2.replace('ener_spec_', '').replace('_', '\_').replace('.dat', ''))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/ener_spec', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/ener_spec', _model,
                                     _file1.replace(".dat", "")+"_comparison.png"),
                        bbox_inches='tight', dpi=300)
        plt.show()


if __name__ == "__main__":

    model = 'heisenberg'
    file1 = 'ener_spec_heisenberg_L_8_obc_J_1_1_1_W_0.5.dat'
    file2 = 'ener_spec_heisenberg_L_8_obc_J_1_1_1_W_8.dat'

    plot_ener_spec(model, file1, file2, _save=True)
