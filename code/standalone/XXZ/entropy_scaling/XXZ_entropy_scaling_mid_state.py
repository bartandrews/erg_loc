import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib import ticker
import time
from joblib import delayed, Parallel
import matplotlib.colors as colors

from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')

# Hamiltonian parameters
J_x_0 = 1
J_z_0 = 1
W_1, W_2 = 0.5, 8
L_list = [8, 10, 12, 14]
# iteration parameters
numb_itr = 20  # 20000 for L=8,10 or 1000 for L=12,14
numb_jobs = 4  # number of spawned processes used for parallelization


# compute H
def realization(itr, W_val):

    print(f"Iteration {itr + 1} of {numb_itr}")

    S_over_L = []

    for L in L_list:

        print(f"L = {L} of {max(L_list)}")
        t_0 = time.time()

        basis = spin_basis_1d(L)

        J_x = [[J_x_0, i, i+1] for i in range(L-1)]
        J_y = [[J_x_0, i, i+1] for i in range(L-1)]
        J_z = [[J_z_0, i, i+1] for i in range(L-1)]
        h_z = [[np.random.uniform(-W_val, W_val), i] for i in range(L)]
        static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z]]
        dynamic = []
        H = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)

        E = H.eigvalsh()
        E_mid, psi_mid = H.eigsh(k=1, sigma=E[2**L//2], maxiter=1E4)

        # E, psi = H.eigh()
        # E_mid = np.sort(E)[len(E)//2]
        # psi_mid = psi[:, np.argsort(E)[len(E)//2]]

        S_over_L.append(basis.ent_entropy(psi_mid, sub_sys_A=range(basis.L//2))["Sent_A"]/L)
        print(f"Total time (seconds): {int(time.time() - t_0)}")

    return S_over_L


fig = plt.figure(figsize=(10, 5))
ax0 = plt.subplot(111)
ax0.set_title(f"$H=\sum_i \sigma^x_{{i}} \sigma^x_{{i+1}} +\sigma^y_i \sigma^y_{{i+1}} + J \sum_i \sigma^z_i \sigma^z_{{i+1}} + \sum_i h_i \sigma^z_i$ with $h_i\in[-W,W]$ and $J={J_z_0}$ (mid state, {numb_itr} disorders)")

W_list = [0, 2.5, 5]
colors = ['tab:red', 'tab:green', 'tab:purple']
for i, W in enumerate(W_list):
    S_over_L = np.asarray(Parallel(n_jobs=numb_jobs)(delayed(realization)(i, W) for i in range(numb_itr)))
    ax0.plot(L_list, np.mean(S_over_L, axis=0), '.-', c=colors[i], label=f'${W}$')

ax0.set_xlabel('$L$')
ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
# ax0.set_xscale('log', basex=2)
ax0.set_ylabel('$S/L$')
ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax0.set_yscale('log')
ax0.legend(title="$W$", ncol=6)

plt.savefig(f"/home/bart/Documents/papers/MBF/XXZ/entropy_scaling/XXZ_entropy_scaling_mid_state_J_{J_z_0}.png", bbox_inches='tight', dpi=300)
plt.show()
