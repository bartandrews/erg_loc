import time
import numpy as np
import matplotlib.pyplot as plt

from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d
from quspin.tools.Floquet import Floquet
from joblib import delayed, Parallel

# system parameters
L = 8
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
# iteration parameters
numb_itr = 12  # 20000 for L=8,10 or 1000 for L=12,14
numb_N = 31
numb_jobs = 12  # number of spawned processes used for parallelization (-1)


def drive(t, Omega):
    return np.sign(np.cos(Omega*t))


def realization(itr):

    print(f"Iteration {itr+1} of {numb_itr}")

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

    # --- energy absorbed under driving
    E_Tinf = H.trace(time=0) / (2**L)
    E, phi = H.eigh(time=0)

    E_0 = np.min(E)
    phi_0 = phi[:, np.argmin(E)]

    t_list = np.array([0.0, T/4.0, 3.0*T/4.0]) + np.finfo(float).eps
    dt_list = np.array([T/4.0, T/2.0, T/4.0])
    Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, UF=True, thetaF=False)
    # Floq = Floquet({'H': H, 'T': T}, UF=True, thetaF=True)
    UF = Floq.UF

    # thetaF = Floq.thetaF
    # print("thetaF = ", [np.angle(i) for i in thetaF])
    # print("UF = ", UF)
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # mat = ax.matshow(np.abs(UF), cmap='Greys')
    # ax.set_xlabel("$i$")
    # ax.set_ylabel("$j$")
    # ax.set_title(f"$T={T}$, $L={L}$")
    # cbar = plt.colorbar(mat)
    # cbar.set_label("$|U_{\mathrm{F}, i, j}|$")
    # plt.show()
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot([np.angle(i)/np.pi for i in thetaF], '.')
    # ax.set_xlabel("$i$")
    # ax.set_ylabel("$\\theta_{\mathrm{F},i}/\pi$")
    # ax.set_title(f"$T={T}$, $L={L}$")
    # plt.show()
    #
    # 1/0

    phi_N = phi_0
    Q_N = np.zeros(numb_N)
    for n in range(numb_N):

        if n > 0:
            phi_N = UF.dot(phi_N)

        Q_N[n] = (np.real(H.matrix_ele(phi_N, phi_N, time=0)) - E_0) / (E_Tinf - E_0)

    return Q_N


if __name__ == '__main__':

    t_0 = time.time()
    Q_N = np.asarray(Parallel(n_jobs=numb_jobs)(delayed(realization)(i) for i in range(numb_itr)))
    print(f"Total time (seconds): {int(time.time() - t_0)}")
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
    ax.set_title(f"$T={T}$, $L={L}$, averaged over {numb_itr} disorders")
    ax.axhline(1, c='r')
    plt.show()
