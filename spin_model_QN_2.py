import time
import numpy as np
import matplotlib.pyplot as plt

from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d
from quspin.tools.Floquet import Floquet
from joblib import delayed, Parallel

# system parameters
L = 8
N_sites = 2**L
basis = spin_basis_1d(L, a=2)
# Hamiltonian parameters
g_1 = 1
J_1 = -1
g_2 = 1
J_2 = -1
W = 100
# driving parameters
T = 10.0
T_d = T/4.0
Omega = 2*np.pi/T
drive_args = []
t = T/2  # 0 < t < T/2
# iteration parameters
numb_itr = 1  # 20000 for L=8,10 or 1000 for L=12,14
numb_N = 31
numb_jobs = 1  # number of spawned processes used for parallelization (-1)


def drive_H_d(t):

    if -T_d <= t < 0:
        return 1
    elif T <= t < T+T_d:
        return -1
    else:
        return 0


def drive_H_1(t):

    if 0 <= t < T/4 or 3*T/4 <= t < T:
        return 1
    else:
        return 0


def drive_H_2(t):

    if T/4 <= t < 3*T/4:
        return 1
    else:
        return 0


def realization(itr):

    print(f"Iteration {itr+1} of {numb_itr}")

    # compute H_mag
    h_z = [[1, i] for i in range(L)]
    static = [["z", h_z]]
    dynamic = []
    H_mag = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)

    # compute H
    G_1 = [[g_1, i] for i in range(L)]
    J_1_x = [[J_1, i, (i+1) % L] for i in range(0, L, 2)]
    J_1_y = [[J_1, i, (i+1) % L] for i in range(0, L, 2)]
    J_1_z = [[J_1, i, (i+1) % L] for i in range(0, L, 2)]
    G_2 = [[g_2, i] for i in range(L)]
    J_2_x = [[J_2, i, (i-1) % L] for i in range(0, L, 2)]
    J_2_y = [[J_2, i, (i-1) % L] for i in range(0, L, 2)]
    J_2_z = [[J_2, i, (i-1) % L] for i in range(0, L, 2)]
    h_z = [[np.random.uniform(-W, W), i] for i in range(L)]
    static = []
    dynamic = [["z", h_z, drive_H_d, drive_args],
               ["I", G_1, drive_H_1, drive_args], ["xx", J_1_x, drive_H_1, drive_args],
               ["yy", J_1_y, drive_H_1, drive_args], ["zz", J_1_z, drive_H_1, drive_args],
               ["I", G_2, drive_H_2, drive_args], ["xx", J_2_x, drive_H_2, drive_args],
               ["yy", J_2_y, drive_H_2, drive_args], ["zz", J_2_z, drive_H_2, drive_args]]
    H = Omega*hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)

    # --- magentization
    if 0 <= t < T/4:
        t_list = np.array([-T_d, 0.0]) + np.finfo(float).eps
        dt_list = np.array([T_d, t])
    elif T/4 <= t < 3*T/4:
        t_list = np.array([-T_d, 0.0, T/4.0]) + np.finfo(float).eps
        dt_list = np.array([T_d, T/4.0, t-T/4.0])
    elif 3*T/4 <= t < T:
        t_list = np.array([-T_d, 0.0, T/4.0, 3.0*T/4.0]) + np.finfo(float).eps
        dt_list = np.array([T_d, T/4.0, T/2.0, t-3.0*T/4.0])
    elif T <= t <= T + T_d:
        t_list = np.array([-T_d, 0.0, T/4.0, 3.0*T/4.0, T]) + np.finfo(float).eps
        dt_list = np.array([T_d, T/4.0, T/2.0, T/4.0, t-T])
    else:
        raise ValueError("Invalid value of t")
    Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, UF=True, thetaF=True)
    # Floq = Floquet({'H': H, 'T': T}, UF=True, thetaF=True)
    UF = Floq.UF

    thetaF = Floq.thetaF
    print("thetaF = ", [np.angle(i) for i in thetaF])
    print("UF = ", UF)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    mat = ax.matshow(np.abs(UF), cmap='Greys')
    ax.set_xlabel("$i$")
    ax.set_ylabel("$j$")
    ax.set_title(f"$t={t}$, $T={T}$, $W={W}$, $L={L}$")
    cbar = plt.colorbar(mat)
    cbar.set_label("$|U_{\mathrm{F}, i, j}|$")
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot([np.angle(i)/np.pi for i in thetaF], '.')
    ax.set_xlabel("$i$")
    ax.set_ylabel("$\\theta_{\mathrm{F},i}/\pi$")
    ax.set_title(f"$t={t}$, $T={T}$, $W={W}$, $L={L}$")
    plt.show()

    1/0

    # phi_0 = np.zeros(N_sites)
    # phi_0[0] = 1
    # phi_N = phi_0
    # mag = np.zeros(numb_N)
    # for n in range(numb_N):
    #
    #     if n > 0:
    #         phi_N = UF.dot(phi_N)
    #
    #     mag[n] = np.real(H_mag.matrix_ele(phi_N, phi_N)) / L

    return mag


if __name__ == '__main__':

    t_0 = time.time()
    mag = np.asarray(Parallel(n_jobs=numb_jobs)(delayed(realization)(i) for i in range(numb_itr)))
    print(f"Total time (seconds): {int(time.time() - t_0)}")
    print(f"Average time per iteration (seconds): {(time.time()-t_0)/numb_itr:.1f}")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    N = np.arange(numb_N)
    # plot the iterations
    for itr in range(numb_itr):
        ax.plot(N, mag[itr], '-', lw=0.1)
    # plot the average
    ax.plot(N, np.mean(mag, axis=0), '-', marker='x', c='k', lw=1)
    ax.set_xlabel("$N$")
    ax.set_ylabel("$\langle m_z \\rangle$")
    plt.rcParams['axes.titlepad'] = 14
    ax.set_title(f"$t={t}$, $T={T}$, $W={W}$, $L={L}$, {numb_itr} disorders")
    ax.axhline(1, c='r')
    plt.show()
