import numpy as np
import matplotlib.pyplot as plt

from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d
from quspin.tools.Floquet import Floquet

# system parameters
L = 10
basis = spin_basis_1d(L, a=2)
# Hamiltonian parameters
g_1 = 1
J_1 = -1
g_2 = 1
J_2 = -1
# driving parameters
T = 10
Omega = 2*np.pi/T
drive_args = [Omega]


def drive(t, Omega):
    return np.sign(np.cos(Omega*t))


# compute H
G_1 = [[g_1, i] for i in range(L)]
J_1_x = [[J_1, i, (i+1) % L] for i in range(0, L, 2)]
J_1_y = [[J_1, i, (i+1) % L] for i in range(0, L, 2)]
J_1_z = [[J_1, i, (i+1) % L] for i in range(0, L, 2)]
pos_G_2 = [[g_2, i] for i in range(L)]
pos_J_2_x = [[J_2, i, (i-1) % L] for i in range(0, L, 2)]
pos_J_2_y = [[J_2, i, (i-1) % L] for i in range(0, L, 2)]
pos_J_2_z = [[J_2, i, (i-1) % L] for i in range(0, L, 2)]
neg_G_2 = [[-g_2, i] for i in range(L)]
neg_J_2_x = [[-J_2, i, (i-1) % L] for i in range(0, L, 2)]
neg_J_2_y = [[-J_2, i, (i-1) % L] for i in range(0, L, 2)]
neg_J_2_z = [[-J_2, i, (i-1) % L] for i in range(0, L, 2)]
static = [["I", G_1], ["xx", J_1_x], ["yy", J_1_y], ["zz", J_1_z],
          ["I", pos_G_2], ["xx", pos_J_2_x], ["yy", pos_J_2_y], ["zz", pos_J_2_z]]
dynamic = [["I", G_1, drive, drive_args], ["xx", J_1_x, drive, drive_args],
           ["yy", J_1_y, drive, drive_args], ["zz", J_1_z, drive, drive_args],
           ["I", neg_G_2, drive, drive_args], ["xx", neg_J_2_x, drive, drive_args],
           ["yy", neg_J_2_y, drive, drive_args], ["zz", neg_J_2_z, drive, drive_args]]
H = 0.5*Omega*hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)

E_Tinf = H.trace(time=0) / (2**L)
E, alpha = H.eigh(time=0)

t_list = np.array([0.0, T/4.0, 3.0*T/4.0]) + np.finfo(float).eps
dt_list = np.array([T/4.0, T/2.0, T/4.0])
Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, VF=True)
psi = Floq.VF

A2 = np.zeros((len(psi), len(alpha)))
for i_idx in range(len(psi)):
    for alpha_idx in range(len(alpha)):
        A2[alpha_idx, i_idx] = np.abs(np.dot(psi[:, i_idx], alpha[:, alpha_idx]))**2

fig = plt.figure()
ax = fig.add_subplot(111)
# for i_idx in range(len(psi)):
#     ax.plot(range(len(alpha)), A2[:, i_idx], '.', label=f"$i={i_idx}$")
ax.plot(range(len(alpha)), A2[:, 0], '.')
ax.set_xlabel("$\\alpha$")
ax.set_ylabel("$|A_{\\alpha, 0}|^2$")
ax.set_yscale('log')
ax.set_ylim([10e-15, 1])
ax.set_title(f"$L={L}$")
ax.axhline(0.0001, c='r')
plt.show()
