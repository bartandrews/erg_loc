import time
import sys
from scipy import signal
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
W = 0.5  # W=0.5 for ergodic or W=8 for MBL
h_0 = 2
# driving parameters
T_1 = 1.5
T_0 = 7
T = T_1 + T_0
drive_args = []
# iteration parameters
numb_itr = 12  # 20000 for L=8,10 or 1000 for L=12,14
numb_N = 31
numb_jobs = -1  # number of spawned processes used for parallelization


def drive_start_V(t):

    # bring t into the range 0 <= t < T
    if t > 0:
        t = t % T
    elif t < 0:
        t = T - abs(t) % T
    else:
        return 0

    # assign a sign to t value
    if 0 < t < T_1:
        return -1  # yields V
    elif T_1 < t < T:
        return +1  # yields H_0
    else:
        return 0  # intermediate values of step function

    # return -signal.square(2*np.pi*t/T, duty=T_1/T)  # period T, step at T_1, V then H_0
    # return np.cos(2*np.pi*t/T) + np.cos(2*np.pi*t/T_1)


def drive_start_H0(t):

    # bring t into the range 0 < t < T
    if t > 0:
        t = t % T
    elif t < 0:
        t = T - abs(t) % T
    else:
        return 0

    # assign a sign to t value
    if 0 < t < T_0:
        return +1  # yields H_0
    elif T_0 < t < T:
        return -1  # yields V
    else:
        return 0


def drive_mid_H0(t):

    # bring t into the range 0 < t < T
    if t >= 0:  # = sign here?
        t = t % T
    elif t < 0:
        t = T - abs(t) % T
    else:
        return 0

    # assign a sign to t value
    if 0 <= t < T_0/2:  # = sign here?
        return +1  # yields H_0
    elif T_0/2 < t < T_0/2 + T_1:
        return -1  # yields V
    elif T_0/2 + T_1 < t < T:
        return +1  # yields H_0
    else:
        return 0


# # plot drive(t)
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ts = np.linspace(0, 2*T, 171)
# ax.plot(ts, [drive(i) for i in ts], '.-', lw=1)
# ax.set_xlabel("t")
# ax.set_xticks(np.arange(int(np.min(ts))-3, int(np.max(ts))+4, 2))
# ax.set_ylabel("drive(t)")
# ax.set_title(f"$T_0={T_0}$, $T_1={T_1}$")
# plt.show()
# 1/0


# compute H
J_x = [[J_x_0, i, i+1] for i in range(L-1)]
J_y = [[J_x_0, i, i+1] for i in range(L-1)]
J_z = [[J_z_0, i, i+1] for i in range(L-1)]
h_z = [[np.random.uniform(-W, W), i] for i in range(L)]
pos_h = [[+h_0, L//2]]
neg_h = [[-h_0, L//2]]

# static_init = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z]]
# dynamic_init = []
# H_init = hamiltonian(static_init, dynamic_init, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)

static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z], ["z", pos_h]]
dynamic_start_V = [["xx", J_x, drive_start_V, drive_args], ["yy", J_y, drive_start_V, drive_args], ["zz", J_z, drive_start_V, drive_args],
                   ["z", h_z, drive_start_V, drive_args], ["z", neg_h, drive_start_V, drive_args]]
dynamic_start_H0 = [["xx", J_x, drive_start_H0, drive_args], ["yy", J_y, drive_start_H0, drive_args], ["zz", J_z, drive_start_H0, drive_args],
                    ["z", h_z, drive_start_H0, drive_args], ["z", neg_h, drive_start_H0, drive_args]]
dynamic_mid_H0 = [["xx", J_x, drive_mid_H0, drive_args], ["yy", J_y, drive_mid_H0, drive_args], ["zz", J_z, drive_mid_H0, drive_args],
                  ["z", h_z, drive_mid_H0, drive_args], ["z", neg_h, drive_mid_H0, drive_args]]

H_start_V = 0.5 * hamiltonian(static, dynamic_start_V, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)
H_start_H0 = 0.5 * hamiltonian(static, dynamic_start_H0, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)
H_mid_H0 = 0.5 * hamiltonian(static, dynamic_mid_H0, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)

Floq_cont_start_V = Floquet({'H': H_start_V, 'T': T})
Floq_dis_start_V = Floquet({'H': H_start_V,
                            't_list': np.array([0.0, T_1]) + np.finfo(float).eps,
                            'dt_list': np.array([T_1, T_0])})

Floq_cont_start_H0 = Floquet({'H': H_start_H0, 'T': T})
Floq_dis_start_H0 = Floquet({'H': H_start_H0,
                            't_list': np.array([0.0, T_0]) + np.finfo(float).eps,
                            'dt_list': np.array([T_0, T_1])})

Floq_cont_mid_H0 = Floquet({'H': H_mid_H0, 'T': T})
Floq_dis_mid_H0 = Floquet({'H': H_mid_H0,
                            't_list': np.array([0.0, T_0/2.0, T_0/2.0+T_1]) + np.finfo(float).eps,
                            'dt_list': np.array([T_0/2.0, T_1, T_0/2.0])})

# plot drive(t)
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(Floq_cont_start_V.EF, '.-', marker='x', label='Floq_cont_start_V')
ax.plot(Floq_cont_start_H0.EF, '.-', marker='^', label='Floq_cont_start_H0')
ax.plot(Floq_cont_mid_H0.EF, '.-', marker='v', label='Floq_cont_mid_H0')

ax.plot(Floq_dis_start_V.EF, '.-', marker='_', label='Floq_dis_start_V')
ax.plot(Floq_dis_start_H0.EF, '.-', marker='s', label='Floq_dis_start_H0')
ax.plot(Floq_dis_mid_H0.EF, '.-', marker='o', label='Floq_dis_mid_H0')

ax.legend()
ax.set_xlabel("$i$")
ax.set_ylabel("$E$")
ax.set_title(f"$T_0={T_0}$, $T_1={T_1}$")
plt.show()
