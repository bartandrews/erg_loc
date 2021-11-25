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
h_0 = 2
W_1, W_2 = 0.5, 8

# compute H_1
J_x = [[J_x_0, i, i+1] for i in range(L-1)]
J_y = [[J_x_0, i, i+1] for i in range(L-1)]
J_z = [[J_z_0, i, i+1] for i in range(L-1)]
h_z = [[np.random.uniform(-W_1, W_1), i] for i in range(L)]
static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z]]
dynamic = []
H_1 = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False, check_pcon=False)
# eigenenergies of H_1
E_1 = H_1.eigvalsh()

r_1 = []
for i in range(1, len(E_1)-1):
    r_1.append((E_1[i+1]-E_1[i])/(E_1[i]-E_1[i-1]))

# compute H_2
J_x = [[J_x_0, i, i+1] for i in range(L-1)]
J_y = [[J_x_0, i, i+1] for i in range(L-1)]
J_z = [[J_z_0, i, i+1] for i in range(L-1)]
h_z = [[np.random.uniform(-W_2, W_2), i] for i in range(L)]
static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z]]
dynamic = []
H_2 = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False, check_pcon=False)
# eigenenergies of H_2
E_2 = H_2.eigvalsh()

r_2 = []
for i in range(1, len(E_2)-1):
    r_2.append((E_2[i+1]-E_2[i])/(E_2[i]-E_2[i-1]))


fig = plt.figure(figsize=(10, 5))
fig.suptitle("Energy Spacing Statistics, $P(r)$ with $r_n=(E_{n+1}-E_n)/(E_n-E_{n-1})$")
gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0)
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1], sharey=ax0)

# plot ergodic
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
ax0.set_title(f"Ergodic ($W={W_1}$, $L={L}$)")

# plot MBL
ax1.yaxis.set_visible(False)
ax1.hist(r_2, bins=np.arange(0, 5.1, 0.1), density=True)
ax1.plot(r_vals, Poisson, c='r', label='Poisson')
ax1.plot(r_vals, GOE, c='g', label='GOE')
ax1.legend(loc='upper right')
ax1.set_xlabel("$r$")
ax1.set_xlim([0, 5])
ax1.set_title(f"MBL ($W={W_2}$, $L={L}$)")

# plt.savefig(f"/home/bart/Documents/papers/MBF/XXZ/energy_spacing_statistics/XXZ_energy_spacing_statistics_hist_J_{J_z_0}.png", bbox_inches='tight', dpi=300)
plt.show()
