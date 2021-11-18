import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib import ticker
import time
from joblib import delayed, Parallel
from numpy.random import default_rng

from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d
from quspin.tools.Floquet import Floquet

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')

# system parameters
L = 10
basis = spin_basis_1d(L, pauli=False)

# compute H_1
J_x = [[1, i, (i+1) % L] for i in range(L)]
J_y = [[1, i, (i+1) % L] for i in range(L)]
J_z = [[1, i, (i+1) % L] for i in range(L)]
rng = default_rng()
# rnd_number = rng.normal(loc=0, scale=10)
# print(rnd_number)
h_z = [[rng.normal(loc=0, scale=100), i] for i in range(L)]
static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z]]
dynamic = []
H_1 = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)
# eigenbasis of H_1
E_1, alpha_1 = H_1.eigh(time=0)

r_1 = []
for i in range(1, len(E_1)-1):
    r_1.append((E_1[i+1]-E_1[i])/(E_1[i]-E_1[i-1]))

# # compute H_2
# J_x = [[J_x_0, i, i+1] for i in range(L-1)]
# J_y = [[J_x_0, i, i+1] for i in range(L-1)]
# J_z = [[J_z_0, i, i+1] for i in range(L-1)]
# h_z = [[8, i] for i in range(L)]
# static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z]]
# dynamic = []
# H_2 = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)
# # eigenbasis of H_1
# E_2, alpha_2 = H_2.eigh(time=0)
#
# r_2 = []
# for i in range(1, len(E_2)-1):
#     r_2.append((E_2[i+1]-E_2[i])/(E_2[i]-E_2[i-1]))


fig = plt.figure(figsize=(10, 5))
fig.suptitle("Energy Spacing Statistics, $P(r)$ with $r_n=(E_{n+1}-E_n)/(E_n-E_{n-1})$")
gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0)
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1], sharey=ax0)

# plot ergodic
print(f"len(r_1)={len(r_1)}")
ax0.hist(r_1, bins=np.arange(0, 5.1, 0.1), density=True)
r_vals = np.arange(0, 100, 0.1)
Poisson = [1/((1+r)**2) for r in r_vals]
GOE = [(27/8)*((r+r**2)/(1+r+r**2)**(5/2)) for r in r_vals]
ax0.plot(r_vals, Poisson, c='r', label='Poisson')
ax0.plot(r_vals, GOE, c='g', label='GOE')
ax0.legend(loc='upper right')
ax0.set_xlabel("$r$")
ax0.set_xlim([0, 5])
ax0.set_ylabel("$P(r)$")
ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax0.set_title(f"Ergodic ($L={L}$)")

# # plot MBL
ax1.yaxis.set_visible(False)
# print(f"len(r_2)={len(r_2)}")
# ax1.hist(r_2, bins=np.arange(0, 5.1, 0.1), density=True)
# ax1.plot(r_vals, Poisson, c='r', label='Poisson')
# ax1.plot(r_vals, GOE, c='g', label='GOE')
# ax1.legend(loc='upper right')
# ax1.set_xlabel("$r$")
# ax1.set_xlim([0, 5])
# ax1.set_title(f"MBL ($W={W_2}$, $L={L}$)")

# plt.savefig("/home/bart/Documents/papers/MBF/Ponte2015/energy_spacing_statistics.png", bbox_inches='tight', dpi=300)
plt.show()
