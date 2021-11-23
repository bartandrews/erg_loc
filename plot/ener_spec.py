import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_energy_spectrum(_model, _file1, _file2=None, _save=False):

    with open(os.path.join('data/energy_spectrum', _model, _file1), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        J_z_1, E = [], []
        for i, row in enumerate(plots):
            J_z_1.append(float(row[0]))
            E.append(float(row[1]))


    fig = plt.figure(figsize=(5, 5))
    ax0 = plt.subplot(111)

    ax0.plot(J_z_1, E, '.', marker='_', c='k', lw=1)
    ax0.set_xlabel('$J$')
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel('$E$')
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))

    if _save:
        plt.savefig(os.path.join('figures/energy_spectrum', _model, 'heisenberg.png'), bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":

    model = 'heisenberg'
    file1 = 'heisenberg.dat'

    plot_energy_spectrum(model, file1, _save=True)
