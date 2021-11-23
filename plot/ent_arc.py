# --- python imports
import matplotlib.pyplot as plt
import csv
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import os

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')


def plot_ent_arc(_model, _params1, _params2=None, _save=False):

    proj_root = '/home/bart/PycharmProjects/erg_loc'

    ener_file1 = os.path.join(proj_root, 'data/ener_spec', _model, f"ener_spec_{_model}_{_params1}.dat")
    ent_file1 = os.path.join(proj_root, 'data/ent_spec', _model, f"ent_spec_{_model}_{_params1}.dat")

    if os.path.exists(ener_file1) is False or os.path.exists(ent_file1) is False:
        raise ValueError("ener_file1 or ent_file1 do not exist")

    with open(ener_file1, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        E_1 = []
        for i, row in enumerate(plots):
            E_1.append(float(row[0]))
    with open(ent_file1, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        S_1 = []
        for i, row in enumerate(plots):
            S_1.append(float(row[0]))

    if _params2 is None:

        plt.figure(figsize=(5, 5))
        ax0 = plt.subplot(111)
        ax0.plot(E_1, S_1, '.', c='k', lw=1)
        ax0.set_xlabel('$E$')
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel('$S$')
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(f"{_model}_{_params1}".replace('_', '\_'))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/ent_arc', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/ent_arc', _model, f"ent_arc_{_model}_{_params1}.png"),
                        bbox_inches='tight', dpi=300)
        plt.show()

    else:

        ener_file2 = os.path.join(proj_root, 'data/ener_spec', _model, f"ener_spec_{_model}_{_params2}.dat")
        ent_file2 = os.path.join(proj_root, 'data/ent_spec', _model, f"ent_spec_{_model}_{_params2}.dat")

        if os.path.exists(ener_file2) is False or os.path.exists(ent_file2) is False:
            raise ValueError("ener_file2 or ent_file2 do not exist")

        with open(ener_file2, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            E_2 = []
            for i, row in enumerate(plots):
                E_2.append(float(row[0]))
        with open(ent_file2, 'r') as csvfile:
            plots = csv.reader(csvfile, delimiter='\t')
            S_2 = []
            for i, row in enumerate(plots):
                S_2.append(float(row[0]))

        plt.figure(figsize=(10, 5))
        gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0)
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1], sharey=ax0)

        ax0.plot(E_1, S_1, '.', c='k', lw=1)
        ax0.set_xlabel('$E$')
        ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_ylabel('$S$')
        ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax0.set_title(f"{_model}_{_params1}".replace('_', '\_'))

        ax1.yaxis.set_visible(False)
        ax1.plot(E_2, S_2, '.', c='k', lw=1)
        ax1.set_xlabel('$E$')
        ax1.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
        ax1.set_title(f"{_model}_{_params2}".replace('_', '\_'))

        if _save:
            os.makedirs(os.path.join(proj_root, 'figures/ent_arc', _model), exist_ok=True)
            plt.savefig(os.path.join(proj_root, 'figures/ent_arc', _model,
                                     f"ent_arc_{_model}_{_params1}_comparison.png"), bbox_inches='tight', dpi=300)
        plt.show()


if __name__ == "__main__":

    model = 'heisenberg'
    params1 = 'L_12_obc_J_1_1_1_W_0.5'
    params2 = 'L_12_obc_J_1_1_1_W_8'

    plot_ent_arc(model, params1, params2, _save=True)
