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

# system parameters
L = 8
basis = spin_basis_1d(L)
# Hamiltonian parameters
J_x_0 = 1
J_z_0 = 1
W = 3
h_0 = 2
# evolution parameters
t_initial_idx = 1
t_final_idx = 3
t_list = np.zeros(t_final_idx-t_initial_idx+1)
for i, t_idx in enumerate(range(t_initial_idx, t_final_idx+1)):
    t_list[i] = 10**t_idx
# iteration parameters
numb_itr = 1  # 20000 for L=8,10 or 1000 for L=12,14
numb_jobs = -1  # number of spawned processes used for parallelization


# compute H
def realization(itr, W_val):

    print(f"Iteration {itr + 1} of {numb_itr}")

    S = []

    J_x = [[J_x_0, i, i+1] for i in range(L-1)]
    J_y = [[J_x_0, i, i+1] for i in range(L-1)]
    J_z = [[J_z_0, i, i+1] for i in range(L-1)]
    h_z = [[np.random.uniform(-W_val, W_val), i] for i in range(L)]
    static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z]]
    dynamic = []
    H = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False, check_pcon=False)

    E, psi = H.eigh()

    E_mid = np.sort(E)[len(E) // 2]
    psi_mid = psi[:, np.argsort(E)[len(E) // 2]]

    psi = H.evolve(psi_mid, 0.0, t_list)

    for i in range(len(t_list)):
        S.append(basis.ent_entropy(psi[:, i], sub_sys_A=range(basis.L//2))["Sent_A"])

    return S


def Page_value(L_val):

    # value = 0

    n = 2**(L_val//2)
    m = n

    print(n, m)

    value = np.log(m) - m / (2 * n)
    # for k in range(n+1, m*n+1):
    #     value += 1/k - (m-1)/(2*n)

    return value


fig = plt.figure(figsize=(10, 5))
ax0 = plt.subplot(111)
ax0.set_title(f"$H=\sum_i \sigma^x_{{i}} \sigma^x_{{i+1}} +\sigma^y_i \sigma^y_{{i+1}} + J \sum_i \sigma^z_i \sigma^z_{{i+1}} + \sum_i h_i \sigma^z_i$ with $h_i\in[-W,W]$, $J={J_z_0}$, $L={L}$ (mid state, {numb_itr} disorders)")

t_0 = time.time()
S_1 = np.asarray(Parallel(n_jobs=numb_jobs)(delayed(realization)(i, W) for i in range(numb_itr)))
print(f"Total time (seconds): {int(time.time() - t_0)}")

ax0.plot(t_list, np.mean(S_1, axis=0), '.-', marker='x', c='k', lw=1, label=f"$W={W}$")
# ax0.plot([min(t_list), max(t_list)], [Page_value(L), Page_value(L)], label=f"Page value")
ax0.legend()
ax0.set_xlabel('$t$')
ax0.set_xscale('log')
ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax0.set_ylabel('$S$')
ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))

plt.savefig(f"/home/bart/Documents/papers/MBF/XXZ/entropy_spreading/XXZ_entropy_spreading_mid_state_J_{J_z_0}.png", bbox_inches='tight', dpi=300)
plt.show()
