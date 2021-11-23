import time
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
W = 0.5  # W=0.5 for ergodic and W=8 for MBL
h_0 = 2
# driving parameters
T_1 = 1.5
T_0 = 7
T = T_1 + T_0
T_H_0 = T_1 + 0.5*T_0
drive_args = []
# iteration parameters
numb_itr = 3  # 20000 for L=8,10 or 1000 for L=12,14
numb_N = 31
# initialize results arrays
Q_N = np.zeros((numb_itr, numb_N))


def drive(t):

    # bring t into the range 0 < t < T
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
        return 0


# # plot drive(t)
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ts = np.linspace(-2*T, 2*T, 400)
# ax.plot(ts, [drive(i) for i in ts], '.-', lw=1)
# ax.set_xlabel("t")
# ax.set_xticks(np.arange(int(np.min(ts))-3, int(np.max(ts))+4, 2))
# ax.set_ylabel("drive(t)")
# ax.set_title(f"$T_1={T_1}$, $T_0={T_0}$")
# plt.show()


t_0 = time.time()
#
for itr in range(numb_itr):

    print(f"Iteration {itr+1} of {numb_itr}")

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

    # --- energy absorbed under driving
    E_Tinf = H.trace(time=T_H_0) / (2**L)
    E, phi = H.eigh(time=T_H_0)

    E_0 = np.min(E)
    phi_0 = phi[:, np.argmin(E)]

    t_list = np.array([0.0, T_1]) + np.finfo(float).eps
    dt_list = np.array([T_1, T_0])
    Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, UF=True)
    # Floq = Floquet({'H': H, 'T': T_0+T_1}, UF=True)
    UF = Floq.UF

    phi_N = phi_0
    for n in range(numb_N):

        if n > 0:
            phi_N = UF.dot(phi_N)

        Q_N[itr][n] = (np.real(H.matrix_ele(phi_N.conj().T, phi_N, time=T_H_0)) - E_0) / (E_Tinf - E_0)
#
print(f"Total time (seconds): {int(time.time()-t_0)}")
print(f"Average time per iteration (seconds): {(time.time()-t_0)/numb_itr:.1f}")


fig = plt.figure()
ax = fig.add_subplot(111)
N = np.arange(numb_N)
# plot the iterations
for itr in range(numb_itr):
    ax.plot(N, Q_N[itr], '-', lw=0.1)
# plot the average
ax.plot(N, np.mean(Q_N, axis=0), '-', marker='x', c='k', lw=1)
ax.set_xlabel("$N$")
ax.set_ylabel("$Q_N$")
ax.set_title(f"$W={W}$, $L={L}$, averaged over {numb_itr} disorders")
ax.axhline(1, c='r')
plt.show()
