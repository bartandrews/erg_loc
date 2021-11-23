import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib import ticker
import time
from joblib import delayed, Parallel

from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d
from quspin.tools.Floquet import Floquet
from quspin.tools.measurements import ent_entropy, diag_ensemble  # entropies

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')

# Hamiltonian parameters
J_x_0 = 1
# J_z_0 = 0.2
J_z_0_list = np.arange(0, 1.01, 0.1)
W_1, W_2 = 0.5, 8
h_0 = 2
L = 8
# iteration parameters
numb_itr = 1  # 20000 for L=8,10 or 1000 for L=12,14
numb_jobs = -1  # number of spawned processes used for parallelization

# compute H
def realization(itr, W_val):

    print(f"Iteration {itr + 1} of {numb_itr}")

    E = np.zeros((len(J_z_0_list), 2**L))

    for i, J_z_0 in enumerate(J_z_0_list):

        basis = spin_basis_1d(L)

        J_x = [[J_x_0, i, i+1] for i in range(L-1)]
        J_y = [[J_x_0, i, i+1] for i in range(L-1)]
        J_z = [[J_z_0, i, i+1] for i in range(L-1)]
        h_z = [[np.random.uniform(-W_val, W_val), i] for i in range(L)]
        static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z]]
        dynamic = []
        H = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)

        E[i] = H.eigvalsh()

    return E


fig = plt.figure(figsize=(10, 5))
fig.suptitle(f"$H=\sum_i \sigma^x_{{i}} \sigma^x_{{i+1}} +\sigma^y_i \sigma^y_{{i+1}} + J \sum_i \sigma^z_i \sigma^z_{{i+1}} + \sum_i h_i \sigma^z_i$ with $h_i\in[-W,W]$")
gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0)
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1], sharey=ax0)


t_0 = time.time()
E_1 = np.asarray(Parallel(n_jobs=numb_jobs)(delayed(realization)(i, W_1) for i in range(numb_itr)))
E_2 = np.asarray(Parallel(n_jobs=numb_jobs)(delayed(realization)(i, W_2) for i in range(numb_itr)))
print(f"Total time (seconds): {int(time.time() - t_0)}")

# for itr in range(numb_itr):
#     ax0.plot(L_list, S_av_1[itr], '-', lw=0.1)

J_z_0_list_plotting = [val for val in J_z_0_list for _ in range(2**L)]

ax0.plot(J_z_0_list_plotting, np.mean(E_1, axis=0).flatten(), '.', marker='_', c='k', lw=1)
ax0.set_xlabel('$J$')
ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax0.set_ylabel('$E$')
ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax0.set_title(f'Ergodic ($W={W_1}$, {numb_itr} disorders)')

ax1.yaxis.set_visible(False)
# for itr in range(numb_itr):
#     ax1.plot(L_list, S_av_2[itr], '-', lw=0.1)
ax1.plot(J_z_0_list, np.mean(E_2, axis=0), '.', marker='_', c='k', lw=1)
ax1.set_xlabel('$J$')
ax1.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax1.set_ylabel('$E$')
ax1.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax1.set_title(f'MBL ($W={W_2}$, {numb_itr} disorders)')

# plt.savefig(f"/home/bart/Documents/papers/MBF/XXZ/energy_spectrum/XXZ_energy_spectrum.png", bbox_inches='tight', dpi=300)
plt.show()
