import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter

from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d
from quspin.tools.Floquet import Floquet
from joblib import delayed, Parallel

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')

# system parameters
L = 8
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
# iteration parameters
numb_itr = 100  # 20000 for L=8,10 or 1000 for L=12,14
numb_N = 31
numb_jobs = -1  # number of spawned processes used for parallelization


def drive(t):

    # bring t into the range 0 < t < T
    if t > 0:  # = sign here?
        t = t % T
    elif t < 0:
        t = T - abs(t) % T
    else:
        return 0

    # assign a sign to t value
    if 0 < t < T_0/2:  # = sign here?
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


def realization(itr, W_val):

    print(f"Iteration {itr+1} of {numb_itr}")

    # compute H
    J_x = [[J_x_0, i, i+1] for i in range(L-1)]
    J_y = [[J_x_0, i, i+1] for i in range(L-1)]
    J_z = [[J_z_0, i, i+1] for i in range(L-1)]
    h_z = [[np.random.uniform(-W_val, W_val), i] for i in range(L)]
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

    t_list = np.array([0.0, T_0/2.0, T_0/2.0+T_1]) + np.finfo(float).eps
    dt_list = np.array([T_0/2.0, T_1, T_0/2.0])
    Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, UF=True)
    # Floq = Floquet({'H': H, 'T': T}, UF=True)
    UF = Floq.UF

    phi_N = phi_0
    Q_N = np.zeros(numb_N)
    for n in range(numb_N):

        if n > 0:
            phi_N = UF.dot(phi_N)

        Q_N[n] = (np.real(H.matrix_ele(phi_N, phi_N, time=0)) - E_0) / (E_Tinf - E_0)

    return Q_N


if __name__ == '__main__':

    fig = plt.figure(figsize=(10, 5))
    fig.suptitle("Energy Absorbed, $Q_N=(\\braket{\phi_0|H_0|\phi_0}-E_0)/(E_{T=\infty}-E_0)$")
    gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1], sharey=ax0)

    t_0 = time.time()
    Q_N_1 = np.asarray(Parallel(n_jobs=numb_jobs)(delayed(realization)(i, 0.5) for i in range(numb_itr)))
    Q_N_2 = np.asarray(Parallel(n_jobs=numb_jobs)(delayed(realization)(i, 8) for i in range(numb_itr)))
    print(f"Total time (seconds): {int(time.time() - t_0)}")

    N = np.arange(numb_N)
    # plot the iterations
    for itr in range(numb_itr):
        ax0.plot(N, Q_N_1[itr], '-', lw=0.1)
        ax1.plot(N, Q_N_2[itr], '-', lw=0.1)
    # plot the average (ergodic)
    ax0.plot(N, np.mean(Q_N_1, axis=0), '-', marker='x', c='k', lw=1)
    ax0.set_xlabel("$N$")
    ax0.set_ylabel("$Q_N$")
    ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
    ax0.set_title(f"Ergodic ($W=0.5$, $L={L}$, {numb_itr} disorders)")
    ax0.axhline(1, c='r')
    # plot the average (MBL)
    ax1.plot(N, np.mean(Q_N_2, axis=0), '-', marker='x', c='k', lw=1)
    ax1.set_xlabel("$N$")
    ax1.set_ylabel("$Q_N$")
    ax1.set_title(f"MBL ($W=8$, $L={L}$, {numb_itr} disorders)")
    ax1.axhline(1, c='r')

    ax1.yaxis.set_visible(False)

    plt.savefig("/home/bart/Documents/papers/MBF/Ponte_2015/energy_absorbed.png", bbox_inches='tight', dpi=300)
    plt.show()
