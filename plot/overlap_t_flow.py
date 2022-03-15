# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amsfonts}\usepackage{braket}')


def plot_overlap_t_flow(_model, _file1, _file2=None, _file3=None, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    with open(os.path.join(proj_root, 'data/overlap_t_flow', _model, _file1), 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        t_1, overlap_1 = [], []
        for i, row in enumerate(plots):
            if i != 0:  # skip t=0 point
                t_1.append(float(row[0]))
                overlap_1.append(float(row[1]))

    if _file2 is None and _file3 is None:

        plt.figure(figsize=(5, 5))
        ax0 = plt.subplot(111)
        ax0.plot(t_1, overlap_1, '.-', label="mid state")
        ax0.set_xlabel("$t$")
        # ax0.set_xscale('log')
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel("$|\\braket{\\mathbb{Z}_d|e^{-\mathrm{i}Ht}|\\mathbb{Z}_d}|^2$")
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('overlap_t_flow_', '').replace('_', '\_').replace('.dat', ''))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/overlap_t_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/overlap_t_flow', _model,
                                     _file1.replace(".dat", "")+".png"),
                        bbox_inches='tight', dpi=300)
        plt.show()

    elif _file3 is None:

        with open(os.path.join(proj_root, 'data/overlap_t_flow', _model, _file2), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            t_2, overlap_2 = [], []
            for i, row in enumerate(plots):
                if i != 0:  # skip t=0 point
                    t_2.append(float(row[0]))
                    overlap_2.append(float(row[1]))

        plt.figure(figsize=(10, 5))
        gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0)
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1], sharey=ax0)

        ax0.plot(t_1, overlap_1, '.-', label="mid state")
        ax0.set_xlabel("$t$")
        # ax0.set_xscale('log')
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel("$|\\braket{\\mathbb{Z}_d|e^{-\mathrm{i}Ht}|\\mathbb{Z}_d}|^2$")
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('overlap_t_flow_', '').replace('_', '\_').replace('.dat', ''))

        ax1.yaxis.set_visible(False)
        ax1.plot(t_2, overlap_2, '.-', label="mid state")
        ax1.set_xlabel("$t$")
        # ax1.set_xscale('log')
        ax1.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax1.set_title(_file2.replace('overlap_t_flow_', '').replace('_', '\_').replace('.dat', ''))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/overlap_t_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/overlap_t_flow', _model,
                                     _file1.replace(".dat", "")+"_comparison.png"),
                        bbox_inches='tight', dpi=300)
        plt.show()

    else:

        with open(os.path.join(proj_root, 'data/overlap_t_flow', _model, _file2), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            t_2, overlap_2 = [], []
            for i, row in enumerate(plots):
                if i != 0:  # skip t=0 point
                    t_2.append(float(row[0]))
                    overlap_2.append(float(row[1]))

        with open(os.path.join(proj_root, 'data/overlap_t_flow', _model, _file3), 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            t_3, overlap_3 = [], []
            for i, row in enumerate(plots):
                if i != 0:  # skip t=0 point
                    t_3.append(float(row[0]))
                    overlap_3.append(float(row[1]))

        plt.figure(figsize=(10, 5))
        ax0 = plt.subplot(111)
        ax0.plot(t_1, overlap_1, '-', label="$d=2$")
        ax0.plot(t_2, overlap_2, '-', label="$d=3$")
        ax0.plot(t_3, overlap_3, '-', label="$d=4$")
        ax0.set_xlabel("$t$")
        # ax0.set_xscale('log')
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel("$|\\braket{\\mathbb{Z}_d|e^{-\mathrm{i}Ht}|\\mathbb{Z}_d}|^2$")
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(_file1.replace('overlap_t_flow_', '').replace('_', '\_').replace('.dat', ''))

        ax0.legend()

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/overlap_t_flow', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/overlap_t_flow', _model,
                                     _file1.replace(".dat", "") + "_comparison.png"),
                        bbox_inches='tight', dpi=300)
        plt.show()


if __name__ == "__main__":

    model = 'pxp'
    file1 = 'overlap_t_flow_pxp_L_24_obc_t_0_30_1000_J_1_1_1_W_0.dat.Z2'
    file2 = 'overlap_t_flow_pxp_L_24_obc_t_0_30_1000_J_1_1_1_W_0.dat.Z3'
    file3 = 'overlap_t_flow_pxp_L_24_obc_t_0_30_1000_J_1_1_1_W_0.dat.Z4'

    plot_overlap_t_flow(model, file1, file2, file3, _save=True)
