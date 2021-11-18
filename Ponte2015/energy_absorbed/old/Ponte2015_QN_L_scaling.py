import time
import numpy as np
import matplotlib.pyplot as plt

from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d
from quspin.tools.Floquet import Floquet
from joblib import delayed, Parallel

# system parameters
#L = 8

# Hamiltonian parameters
J_x_0 = 1
J_z_0 = 1
W = 0.5  # W=0.5 for ergodic or W=8 for MBL
h_0 = 2
# driving parameters
T_0 = 7
T_1 = 1.5
T = T_0 + T_1
drive_args = []
# iteration parameters
numb_itr = 12  # 20000 for L=8,10 or 1000 for L=12,14
numb_N = 31
numb_jobs = -1  # number of spawned processes used for parallelization


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


# # plot drive(t)
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ts = np.linspace(-2*T, 2*T, 400)
# ax.plot(ts, [drive(i) for i in ts], '.-', lw=1)
# ax.set_xlabel("t")
# ax.set_xticks(np.arange(int(np.min(ts))-3, int(np.max(ts))+4, 2))
# ax.set_ylabel("drive(t)")
# ax.set_title(f"$T_0={T_0}$, $T_1={T_1}$")
# plt.show()


def realization(itr, L):

    print(f"Iteration {itr+1} of {numb_itr}")

    # compute H
    basis = spin_basis_1d(L)
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
    E_Tinf = H.trace(time=0) / (2**L)
    E, phi = H.eigh(time=0)

    E_0 = np.min(E)
    phi_0 = phi[:, np.argmin(E)]

    t_list = np.array([0.0, T_0/2.0, T_0/2+T_1]) + np.finfo(float).eps
    dt_list = np.array([T_0/2.0, T_1, T_0/2.0])
    Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, UF=True)
    # Floq = Floquet({'H': H, 'T': T_0+T_1}, UF=True)
    UF = Floq.UF

    phi_N = phi_0
    Q_N = np.zeros(numb_N)
    for n in range(numb_N):

        if n > 0:
            phi_N = UF.dot(phi_N)

        Q_N[n] = (np.real(H.matrix_ele(phi_N, phi_N, time=0)) - E_0) / (E_Tinf - E_0)

    return Q_N


if __name__ == '__main__':

    L_list = [6, 8, 10]

    Q_N_array = np.zeros((len(L_list), numb_itr, numb_N))

    for length, L_val in enumerate(L_list):
        t_0 = time.time()
        Q_N_array[length] = np.asarray(Parallel(n_jobs=numb_jobs)(delayed(realization)(i, L_val) for i in range(numb_itr)))
        print(f"Total time (seconds): {int(time.time() - t_0)}")
        print(f"Average time per iteration (seconds): {(time.time()-t_0)/numb_itr:.1f}")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    N = np.arange(numb_N)
    # plot the iterations
    # for itr in range(numb_itr):
    #     ax.plot(N, Q_N[itr], '-', lw=0.1)
    # plot the average
    for length, L_val in enumerate(L_list):
        ax.plot(N, np.mean(Q_N_array[length], axis=0), '-', marker='x', lw=1, label=f"$L={L_val}$")
    ax.set_xlabel("$N$")
    ax.set_ylabel("$Q_N$")
    ax.set_title(f"$W={W}$, averaged over {numb_itr} disorders")
    ax.legend(loc='lower right', ncol=5)
    ax.axhline(1, c='r')
    plt.show()
