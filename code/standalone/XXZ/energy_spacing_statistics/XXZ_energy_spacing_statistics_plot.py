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

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')

# system parameters
L = 8
basis = spin_basis_1d(L, pauli=0, Nup=L//2)
# Hamiltonian parameters
J_x_0 = 1
J_z_0 = 1
W_list = np.arange(0.5, 12.5, 0.5)
# iteration parameters
numb_itr = 100  # 20000 for L=8,10 or 1000 for L=12,14
numb_jobs = -1  # number of spawned processes used for parallelization


def realization(itr):

    print(f"Iteration {itr + 1} of {numb_itr}")

    r_av = []
    for W in W_list:

        # compute H
        J_x = [[J_x_0, i, i+1] for i in range(L-1)]
        J_y = [[J_x_0, i, i+1] for i in range(L-1)]
        J_z = [[J_z_0, i, i+1] for i in range(L-1)]
        h_z = [[np.random.uniform(-W, W), i] for i in range(L)]
        static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z]]
        dynamic = []
        H = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False, check_pcon=False)
        # eigenenergies of H
        E = H.eigvalsh()

        r = []
        for i in range(1, len(E)-1):
            delta_n = E[i]-E[i-1]
            delta_n_plus_1 = E[i+1]-E[i]
            r.append(min(delta_n, delta_n_plus_1)/max(delta_n, delta_n_plus_1))
        r_av.append(np.mean(r))

    return r_av


fig = plt.figure(figsize=(10, 5))
ax0 = plt.subplot(111)
ax0.set_title(f"$H=\sum_i \sigma^x_{{i}} \sigma^x_{{i+1}} +\sigma^y_i \sigma^y_{{i+1}} + J \sum_i \sigma^z_i \sigma^z_{{i+1}} + \sum_i h_i \sigma^z_i$ with $h_i\in[-W,W]$, $J={J_z_0}$, $L={L}$ ({numb_itr} disorders)")

t_0 = time.time()
r_av = np.asarray(Parallel(n_jobs=numb_jobs)(delayed(realization)(i) for i in range(numb_itr)))
print(f"Total time (seconds): {int(time.time() - t_0)}")

# plot ergodic
ax0.plot(W_list, np.mean(r_av, axis=0), '.-', label=f"$L={L}$")
ax0.plot([min(W_list), max(W_list)], [0.39, 0.39], label="Poisson")
ax0.plot([min(W_list), max(W_list)], [0.53, 0.53], label="GOE")
ax0.legend()
ax0.set_xlabel("$W$")
ax0.set_xscale('log', basex=2)
ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax0.set_ylabel("$[r^{(n)}_\\alpha]$")
ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))

# plt.savefig(f"/home/bart/Documents/papers/MBF/XXZ/energy_spacing_statistics/XXZ_energy_spacing_statistics_plot_J_{J_z_0}.png", bbox_inches='tight', dpi=300)
plt.show()
