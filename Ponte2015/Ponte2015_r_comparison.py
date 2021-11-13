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
L = 10
basis = spin_basis_1d(L)
# Hamiltonian parameters
J_x_0 = 1
J_z_0 = 1
# W = 0.5  # W=0.5 for ergodic or W=8 for MBL
h_0 = 2
# driving parameters
T_0 = 7
T_1 = 1.5
T = T_0 + T_1
drive_args = []


def drive(t):

    # bring t into the range 0 < t < T
    if t > 0:
        t = t % T
    elif t < 0:
        t = T - abs(t) % T
    else:
        return 0

    # assign a sign to t value
    if 0 < t < T_0/2:
        return +1  # yields H_0
    elif T_0/2 < t < T_0/2 + T_1:
        return -1  # yields V
    elif T_0/2 + T_1 < t < T:
        return +1  # yields H_0
    else:
        return 0


# compute H_1
J_x = [[J_x_0, i, i+1] for i in range(L-1)]
J_y = [[J_x_0, i, i+1] for i in range(L-1)]
J_z = [[J_z_0, i, i+1] for i in range(L-1)]
h_z = [[-0.1, i] for i in range(L)]
pos_h = [[+h_0, L//2]]
neg_h = [[-h_0, L//2]]
static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z], ["z", pos_h]]
dynamic = [["xx", J_x, drive, drive_args], ["yy", J_y, drive, drive_args], ["zz", J_z, drive, drive_args],
           ["z", h_z, drive, drive_args], ["z", neg_h, drive, drive_args]]
H_1 = 0.5*hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)
# eigenbasis of H_1
E_1, alpha_1 = H_1.eigh(time=0)

r_1 = []
for i in range(1, len(E_1)-1):
    r_1.append((E_1[i+1]-E_1[i])/(E_1[i]-E_1[i-1]))

# compute H_2
J_x = [[J_x_0, i, i+1] for i in range(L-1)]
J_y = [[J_x_0, i, i+1] for i in range(L-1)]
J_z = [[J_z_0, i, i+1] for i in range(L-1)]
h_z = [[8, i] for i in range(L)]
pos_h = [[+h_0, L//2]]
neg_h = [[-h_0, L//2]]
static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z], ["z", pos_h]]
dynamic = [["xx", J_x, drive, drive_args], ["yy", J_y, drive, drive_args], ["zz", J_z, drive, drive_args],
           ["z", h_z, drive, drive_args], ["z", neg_h, drive, drive_args]]
H_2 = 0.5*hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)
# eigenbasis of H_1
E_2, alpha_2 = H_2.eigh(time=0)

r_2 = []
for i in range(1, len(E_2)-1):
    r_2.append((E_2[i+1]-E_2[i])/(E_2[i]-E_2[i-1]))


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
ax0.set_title(f"Ergodic ($W=-0.1$, $L={L}$)")

# plot MBL
ax1.yaxis.set_visible(False)
print(f"len(r_2)={len(r_2)}")
ax1.hist(r_2, bins=np.arange(0, 5.1, 0.1), density=True)
ax1.plot(r_vals, Poisson, c='r', label='Poisson')
ax1.plot(r_vals, GOE, c='g', label='GOE')
ax1.legend(loc='upper right')
ax1.set_xlabel("$r$")
ax1.set_xlim([0, 5])
ax1.set_title(f"MBL ($W=8$, $L={L}$)")

plt.savefig("/home/bart/Documents/papers/MBF/Ponte_2015/energy_spacing_statistics.png", bbox_inches='tight', dpi=300)
plt.show()
