# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import os
from scipy.stats import gaussian_kde
import numpy as np
import re

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{amsfonts}\usepackage{braket}')


def plot_overlap_spec(_model, _params1, _params2=None, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    ener_file1 = os.path.join(proj_root, 'data/ener', _model, f"ener_{_model}_{_params1}.dat")
    overlap_file1 = os.path.join(proj_root, 'data/overlap', _model, f"overlap_{_model}_{_params1}.dat")

    if os.path.exists(ener_file1) is False or os.path.exists(overlap_file1) is False:
        raise ValueError("ener_file1 or overlap_file1 do not exist")

    with open(ener_file1, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        E_1 = []
        for i, row in enumerate(plots):
            E_1.append(float(row[0]))
    with open(overlap_file1, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        overlap_1 = []
        for i, row in enumerate(plots):
            overlap_1.append(float(row[0]))

    if _params2 is None:

        plt.figure(figsize=(5, 5))
        ax0 = plt.subplot(111)

        # convert lists to arrays
        E_1 = np.array(E_1)
        overlap_1 = np.array(overlap_1)
        # calculate the point density
        E_1_overlap_1 = np.vstack([E_1, overlap_1])
        z = gaussian_kde(E_1_overlap_1)(E_1_overlap_1)
        # sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        E_1, overlap_1, z = E_1[idx], overlap_1[idx], z[idx]
        # scatter plot
        ax0.scatter(E_1, np.log10(overlap_1), c=z, s=10, cmap='plasma')
        # ax0.plot(E_1, np.log10(overlap_1), '.')
        # red crosses
        L_pattern = re.search('L_(\d+)', _params1)
        L = int(L_pattern.group(1))
        old_frac = 0
        for frac in np.linspace(0.0415, 1.0415, L // 2 + 1):
            tmp_E_1_list, tmp_overlap_list = [], []
            for i, E_1_val in enumerate(E_1):
                if min(E_1) + old_frac * (max(E_1) - min(E_1)) <= E_1_val < min(E_1) + frac * (max(E_1) - min(E_1)):
                    tmp_E_1_list.append(E_1_val)
                    tmp_overlap_list.append(overlap_1[i])
            ax0.axvline(min(E_1) + frac * (max(E_1) - min(E_1)), linewidth=0.1)
            ax0.scatter(tmp_E_1_list[np.argmax(tmp_overlap_list)], np.log10(np.max(tmp_overlap_list)), c='r', marker='x', s=40)
            old_frac = frac
        ax0.set_xlabel('$E$')
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel('$\log_{10} |\\braket{\mathbb{Z}_2|\psi}|^2$')
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylim([-10, 0])
        ax0.set_title(f"{_model}_{_params1}".replace('_', '\_'))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/overlap_spec', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/overlap_spec', _model, f"overlap_spec_{_model}_{_params1}.png"),
                        bbox_inches='tight', dpi=300)
        plt.show()

    else:

        ener_file2 = os.path.join(proj_root, 'data/ener', _model, f"ener_{_model}_{_params2}.dat")
        overlap_file2 = os.path.join(proj_root, 'data/overlap', _model, f"overlap_{_model}_{_params2}.dat")

        if os.path.exists(ener_file2) is False or os.path.exists(overlap_file2) is False:
            raise ValueError("ener_file2 or overlap_file2 do not exist")

        with open(ener_file2, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            E_2 = []
            for i, row in enumerate(plots):
                E_2.append(float(row[0]))
        with open(overlap_file2, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            overlap_2 = []
            for i, row in enumerate(plots):
                overlap_2.append(float(row[0]))

        plt.figure(figsize=(10, 5))
        gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0)
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1], sharey=ax0)

        # convert lists to arrays
        E_1 = np.array(E_1)
        overlap_1 = np.array(overlap_1)
        # calculate the point density
        E_1_overlap_1 = np.vstack([E_1, overlap_1])
        z = gaussian_kde(E_1_overlap_1)(E_1_overlap_1)
        # sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        E_1, overlap_1, z = E_1[idx], overlap_1[idx], z[idx]
        # scatter plot
        ax0.scatter(E_1, np.log10(overlap_1), c=z, s=10, cmap='plasma')
        # red crosses
        # L_pattern = re.search('L_(\d+)', _params1)
        # L = int(L_pattern.group(1))
        # old_frac = 0
        # for frac in np.linspace(0.0415, 1.0415, L//2 + 1):
        #     tmp_E_1_list, tmp_overlap_list = [], []
        #     for i, E_1_val in enumerate(E_1):
        #         if min(E_1) + old_frac * (max(E_1) - min(E_1)) <= E_1_val < min(E_1) + frac * (max(E_1) - min(E_1)):
        #             tmp_E_1_list.append(E_1_val)
        #             tmp_overlap_list.append(overlap_1[i])
        #     ax0.axvline(min(E_1) + frac * (max(E_1) - min(E_1)), linewidth=0.1)
        #     ax0.scatter(tmp_E_1_list[np.argmax(tmp_overlap_list)], np.log10(np.max(tmp_overlap_list)), c='r',
        #                 marker='x', s=40)
        #     old_frac = frac
        ax0.set_xlabel('$E$')
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel('$|\\braket{\mathbb{Z}_2|\psi}|^2$')
        ax0.set_ylim([-10, 0])
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(f"{_model}_{_params1}".replace('_', '\_'))

        ax1.yaxis.set_visible(False)
        # convert lists to arrays
        E_2 = np.array(E_2)
        overlap_2 = np.array(overlap_2)
        # calculate the point density
        E_2_overlap_2 = np.vstack([E_2, overlap_2])
        z = gaussian_kde(E_2_overlap_2)(E_2_overlap_2)
        # sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        E_2, overlap_2, z = E_2[idx], overlap_2[idx], z[idx]
        # scatter plot
        ax1.scatter(E_2, np.log10(overlap_2), c=z, s=10, cmap='plasma')
        # red crosses
        # L_pattern = re.search('L_(\d+)', _params1)
        # L = int(L_pattern.group(1))
        # old_frac = 0
        # for frac in np.linspace(0.0415, 1.0415, L//2 + 1):
        #     tmp_E_2_list, tmp_overlap_list = [], []
        #     for i, E_2_val in enumerate(E_2):
        #         if min(E_2) + old_frac * (max(E_2) - min(E_2)) <= E_2_val < min(E_2) + frac * (max(E_2) - min(E_2)):
        #             tmp_E_2_list.append(E_2_val)
        #             tmp_overlap_list.append(overlap_2[i])
        #     ax1.axvline(min(E_2) + frac * (max(E_2) - min(E_2)), linewidth=0.1)
        #     ax1.scatter(tmp_E_2_list[np.argmax(tmp_overlap_list)], np.log10(np.max(tmp_overlap_list)), c='r',
        #                 marker='x', s=40)
        #     old_frac = frac
        ax1.set_xlabel('$E$')
        ax1.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax1.set_title(f"{_model}_{_params2}".replace('_', '\_'))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/overlap_spec', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/overlap_spec', _model,
                                     f"overlap_spec_{_model}_{_params1}_comparison.png"), bbox_inches='tight', dpi=300)
        plt.show()


if __name__ == "__main__":

    model = 'pxp'
    params1 = 'L_26_obc_J_1_1_1_W_0'
    params2 = 'L_28_obc_J_1_1_1_W_0'

    plot_overlap_spec(model, params1, params2, _save=True)
