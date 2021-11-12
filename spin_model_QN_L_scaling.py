import time
import numpy as np
import matplotlib.pyplot as plt

from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d
from quspin.tools.Floquet import Floquet
from joblib import delayed, Parallel

# system parameters
# L = 10
# Hamiltonian parameters
g_1 = 1
J_1 = [-1, -1, -1]
g_2 = 1
J_2 = [-1, -1, -1]
# driving parameters
T = 10
Omega = 2*np.pi/T
drive_args = [Omega]
# iteration parameters
numb_itr = 12*1  # 20000 for L=8,10 or 1000 for L=12,14
numb_N = 31
numb_jobs = -1  # number of spawned processes used for parallelization


def drive(t, Omega):
    return np.sign(np.cos(Omega*t))


def realization(itr, L):

    print(f"Iteration {itr+1} of {numb_itr}")

    # compute H
    basis = spin_basis_1d(L)
    G_1 = [[g_1, i] for i in range(L)]
    J_x_1 = [[J_1[0], i, i+2] for i in range(L-2)]
    J_y_1 = [[J_1[1], i, i+2] for i in range(L-2)]
    J_z_1 = [[J_1[2], i, i+2] for i in range(L-2)]
    pos_G_2 = [[g_2, i] for i in range(L)]
    pos_J_x_2 = [[J_2[0], i+1, i+2] for i in range(L-2)]
    pos_J_y_2 = [[J_2[1], i+1, i+2] for i in range(L-2)]
    pos_J_z_2 = [[J_2[2], i+1, i+2] for i in range(L-2)]
    neg_G_2 = [[-g_2, i] for i in range(L)]
    neg_J_x_2 = [[-J_2[0], i+1, i+2] for i in range(L-2)]
    neg_J_y_2 = [[-J_2[1], i+1, i+2] for i in range(L-2)]
    neg_J_z_2 = [[-J_2[2], i+1, i+2] for i in range(L-2)]
    static = [["I", G_1], ["xx", J_x_1], ["yy", J_y_1], ["zz", J_z_1],
              ["I", pos_G_2], ["xx", pos_J_x_2], ["yy", pos_J_y_2], ["zz", pos_J_z_2]]
    dynamic = [["I", G_1, drive, drive_args], ["xx", J_x_1, drive, drive_args], ["yy", J_y_1, drive, drive_args], ["zz", J_z_1, drive, drive_args],
               ["I", neg_G_2, drive, drive_args], ["xx", neg_J_x_2, drive, drive_args], ["yy", neg_J_y_2, drive, drive_args], ["zz", neg_J_z_2, drive, drive_args]]
    H = 0.5*Omega*hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)

    # --- energy absorbed under driving
    E_Tinf = H.trace(time=0) / (2**L)
    E, phi = H.eigh(time=0)

    E_0 = np.min(E)
    phi_0 = phi[:, np.argmin(E)]

    t_list = np.array([0.0, T/4.0, 3*T/4.0]) + np.finfo(float).eps
    dt_list = np.array([T/4.0, T/2.0, T/4.0])
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
    ax.set_title(f"$T={T}$, averaged over {numb_itr} runs")
    ax.legend(loc='lower right', ncol=5)
    ax.axhline(1, c='r')
    plt.show()
