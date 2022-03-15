# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import os
import numpy as np
from statistics import mean
from itertools import groupby

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_spec_func(_model, _file, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    with open(os.path.join(proj_root, 'data/spec_func', _model, _file), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        nrow = len(list(plots))
        csvfile.seek(0)

        data = np.zeros((nrow, 2))
        for i, row in enumerate(plots):
            data[i, 0] = float(row[0])
            data[i, 1] = float(row[1])

    plt.figure(figsize=(10, 5))
    ax0 = plt.subplot(111)

    data = data[data[:, 0].astype(float).argsort()]  # sort the data by first column

    ax0.plot(data[:, 0], data[:, 1], '.', c='k', label="scatter", markersize=0.1)
    ax0.axvline(-2.66 * 2, ls='--', c='k')
    ax0.axvline(-2.66, ls='--', c='k')
    ax0.axvline(2.66, ls='--', c='k')
    ax0.axvline(2.66 * 2, ls='--', c='k')

    for i in range(np.shape(data)[0]):  # round the numbers in first column
        data[i, 0] = round(data[i, 0], 2)

    print("before:", data)

    grouper = groupby(data, key=lambda x: x[0])
    data = np.array([[x, mean(yi[1] for yi in y)] for x, y in grouper])  # take mean of degenerate y-values

    print("after:", data)

    ax0.plot(data[:, 0], data[:, 1], '-', c='r', label="line")

    ax0.set_xlabel("$\omega$")
    ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_ylabel("$f^2(\omega)$")
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_yscale('log')
    ax0.set_title(_file.replace('spec_func_', '').replace('_', '\_').replace('.dat', ''))

    if _save:
        os.makedirs(os.path.join(proj_root, 'figures/spec_func', _model), exist_ok=True)
        plt.savefig(os.path.join(proj_root, 'figures/spec_func', _model,
                                 _file.replace(".dat", "")+".png"),
                    bbox_inches='tight', dpi=300)
    plt.show()


if __name__ == "__main__":

    model = 'pxp'
    file = 'spec_func_pxp_L_26_obc_J_1_1_1_W_0.dat'

    plot_spec_func(model, file, _save=True)
