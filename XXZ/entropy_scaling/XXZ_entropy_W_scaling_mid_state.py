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
J_z_0 = 0.2
h_0 = 2
W_list = np.arange(0.1, 8, 0.01)
# iteration parameters
numb_itr = 1000  # 20000 for L=8,10 or 1000 for L=12,14
numb_jobs = -1  # number of spawned processes used for parallelization


# compute H
def realization(itr, L_val):

    print(f"Iteration {itr + 1} of {numb_itr}")

    basis = spin_basis_1d(L_val)

    S = []

    for W in W_list:

        J_x = [[J_x_0, i, i+1] for i in range(L_val-1)]
        J_y = [[J_x_0, i, i+1] for i in range(L_val-1)]
        J_z = [[J_z_0, i, i+1] for i in range(L_val-1)]
        h_z = [[np.random.uniform(-W, W), i] for i in range(L_val)]
        static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z]]
        dynamic = []
        H = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)

        E, psi = H.eigsh()
        E_mid = np.sort(E)[len(E) // 2]
        psi_mid = psi[:, np.argsort(E)[len(E) // 2]]

        S.append(basis.ent_entropy(psi_mid, sub_sys_A=range(basis.L//2))["Sent_A"])

    return S


fig = plt.figure(figsize=(10, 5))
ax0 = plt.subplot(111)
ax0.set_title(f"$H=\sum_i \sigma^x_{{i}} \sigma^x_{{i+1}} +\sigma^y_i \sigma^y_{{i+1}} + J \sum_i \sigma^x_i \sigma^z_{{i+1}} + \sum_i h_i \sigma^z_i$ with $h_i\in[-W,W]$ and $J={J_z_0}$ (mid state, {numb_itr} disorders)")

for L in [8, 10, 12]:
    S_av = np.asarray(Parallel(n_jobs=numb_jobs)(delayed(realization)(i, L) for i in range(numb_itr)))
    # for itr in range(numb_itr):
    #     ax0.plot(W_list, S_av[itr], '-', lw=0.1)
    ax0.plot(W_list, np.mean(S_av, axis=0), '.-', marker='x', label=f'$L={L}$')

ax0.legend()
ax0.set_xlabel('$W$')
ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax0.set_ylabel('$S$')
ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))

plt.savefig(f"/home/bart/Documents/papers/MBF/XXZ/entropy_scaling/entropy_W_scaling_mid_state_J_{J_z_0}.png", bbox_inches='tight', dpi=300)
plt.show()
