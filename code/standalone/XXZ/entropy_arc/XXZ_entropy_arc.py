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
J_z_0 = 1
W_1, W_2 = 0.5, 8
L = 12


# compute H
def realization(W_val):

    S = []

    basis = spin_basis_1d(L)

    J_x = [[J_x_0, i, i+1] for i in range(L-1)]
    J_y = [[J_x_0, i, i+1] for i in range(L-1)]
    J_z = [[J_z_0, i, i+1] for i in range(L-1)]
    h_z = [[np.random.uniform(-W_val, W_val), i] for i in range(L)]
    # h_z = [[8, i] for i in range(L)]
    static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z]]
    dynamic = []
    H = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)

    E, psi = H.eigh()

    print(E.shape)
    print(psi.shape)

    for i in range(psi.shape[1]):
        S.append(basis.ent_entropy(psi[:, i], sub_sys_A=range(basis.L//2))["Sent_A"])

    return E, S


fig = plt.figure(figsize=(10, 5))
fig.suptitle(f"$H=\sum_i \sigma^x_{{i}} \sigma^x_{{i+1}} +\sigma^y_i \sigma^y_{{i+1}} + J \sum_i \sigma^z_i \sigma^z_{{i+1}} + \sum_i h_i \sigma^z_i$ with $h_i\in[-W,W]$, $J={J_z_0}$, $L={L}$")
gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0)
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1], sharey=ax0)


t_0 = time.time()
E_1, S_1 = realization(W_1)
E_2, S_2 = realization(W_2)
print(E_1)
print(E_2)
print(S_1)
print(S_2)
print(f"Total time (seconds): {int(time.time() - t_0)}")

ax0.plot(E_1, S_1, '.', c='k', lw=1)
ax0.set_xlabel('$E$')
ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax0.set_ylabel('$S$')
ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax0.set_title(f'Ergodic ($W={W_1}$)')

ax1.yaxis.set_visible(False)
ax1.plot(E_2, S_2, '.', c='k', lw=1)
ax1.set_xlabel('$E$')
ax1.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax1.set_ylabel('$S$')
ax1.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax1.set_title(f'MBL ($W={W_2}$)')

# plt.savefig(f"/home/bart/Documents/papers/MBF/XXZ/entropy_arc/XXZ_entropy_arc_J_{J_z_0}.png", bbox_inches='tight', dpi=300)
plt.show()
