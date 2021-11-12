import numpy as np
import matplotlib.pyplot as plt

from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d
from quspin.tools.Floquet import Floquet

# system parameters
L = 8
basis = spin_basis_1d(L)
# Hamiltonian parameters
J_x_0 = 1
J_z_0 = 1
W = 8  # W=0.5 for ergodic or W=8 for MBL
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


# compute H
J_x = [[J_x_0, i, i+1] for i in range(L-1)]
J_y = [[J_x_0, i, i+1] for i in range(L-1)]
J_z = [[J_z_0, i, i+1] for i in range(L-1)]
h_z = [[np.random.uniform(-W, W), i] for i in range(L)]
pos_h = [[+h_0, L//2]]
neg_h = [[-h_0, L//2]]
static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z], ["z", pos_h]]
dynamic = [["xx", J_x, drive, drive_args], ["yy", J_y, drive, drive_args], ["zz", J_z, drive, drive_args],
           ["z", h_z, drive, drive_args], ["z", neg_h, drive, drive_args]]
H = 0.5*hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)

E_Tinf = H.trace(time=0) / (2**L)
E, alpha = H.eigh(time=0)

t_list = np.array([0.0, T_0/2.0, T_0/2.0+T_1]) + np.finfo(float).eps
dt_list = np.array([T_0/2.0, T_1, T_0/2.0])
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
ax.set_title(f"$W={W}$, $L={L}$")
ax.axhline(0.0001, c='r')
plt.show()
