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
h_0 = 2
L_list = [8, 10, 12, 14, 16, 18, 20]
# iteration parameters
numb_itr = 100  # 20000 for L=8,10 or 1000 for L=12,14
numb_jobs = -1  # number of spawned processes used for parallelization


# compute H
def realization(itr, W_val):

    print(f"Iteration {itr + 1} of {numb_itr}")

    S = []

    for L in L_list:

        basis = spin_basis_1d(L)

        J_x = [[J_x_0, i, i+1] for i in range(L-1)]
        J_y = [[J_x_0, i, i+1] for i in range(L-1)]
        J_z = [[J_z_0, i, i+1] for i in range(L-1)]
        h_z = [[np.random.uniform(-W_val, W_val), i] for i in range(L)]
        static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z]]
        dynamic = []
        H = 0.5*hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)

        E, psi = H.eigsh(k=1, which="SA", maxiter=1E4)
        S.append(basis.ent_entropy(psi, sub_sys_A=range(basis.L//2))["Sent_A"])

    return S


fig = plt.figure(figsize=(10, 5))
fig.suptitle("AFM XXX model with random z-field $h_i\in[-W,W]$ (ground state)")
gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0)
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1], sharey=ax0)


t_0 = time.time()
S_av_1 = np.asarray(Parallel(n_jobs=numb_jobs)(delayed(realization)(i, W_1) for i in range(numb_itr)))
S_av_2 = np.asarray(Parallel(n_jobs=numb_jobs)(delayed(realization)(i, W_2) for i in range(numb_itr)))
print(f"Total time (seconds): {int(time.time() - t_0)}")

# for itr in range(numb_itr):
#     ax0.plot(W_list, S_av[itr], '-', lw=0.1)
ax0.plot(L_list, np.mean(S_av_1, axis=0), '.-')
ax0.set_xlabel('$L$')
ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax0.set_ylabel('$S$')
ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax0.set_title(f'Ergodic ($W={W_1}$, {numb_itr} disorders)')

ax1.yaxis.set_visible(False)
ax1.plot(L_list, np.mean(S_av_2, axis=0), '.-')
ax1.set_xlabel('$L$')
ax1.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax1.set_ylabel('$S$')
ax1.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax1.set_title(f'MBL ($W={W_2}$, {numb_itr} disorders)')

plt.savefig(f"/home/bart/Documents/papers/MBF/Ponte_2015/entropy_scaling/entropy_L_scaling_Hamiltonian.png", bbox_inches='tight', dpi=300)
plt.show()
