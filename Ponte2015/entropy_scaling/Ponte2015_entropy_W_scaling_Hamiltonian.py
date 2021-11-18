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
from quspin.tools.measurements import ent_entropy, diag_ensemble  # entropies


plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}\usepackage{braket}')

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
W_list = np.arange(0.1, 8, 0.01)
#
numb_itr = 1000
numb_jobs = 12


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
def realization(itr, L_val):

    print(f"Iteration {itr + 1} of {numb_itr}")

    basis = spin_basis_1d(L_val)

    S = []

    for W in W_list:

        J_x = [[J_x_0, i, i+1] for i in range(L_val-1)]
        J_y = [[J_x_0, i, i+1] for i in range(L_val-1)]
        J_z = [[J_z_0, i, i+1] for i in range(L_val-1)]
        h_z = [[np.random.uniform(-W, W), i] for i in range(L_val)]
        pos_h = [[+h_0, L_val//2]]
        neg_h = [[-h_0, L_val//2]]
        static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z], ["z", pos_h]]
        dynamic = [["xx", J_x, drive, drive_args], ["yy", J_y, drive, drive_args], ["zz", J_z, drive, drive_args],
                   ["z", h_z, drive, drive_args], ["z", neg_h, drive, drive_args]]
        H = 0.5*hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)

        E, psi = H.eigsh(k=1, which="SA", maxiter=1E4)
        S.append(basis.ent_entropy(psi, sub_sys_A=range(basis.L//2))["Sent_A"])

    return S


fig = plt.figure(figsize=(10, 5))
ax0 = plt.subplot(111)

for L in [12]:
    S_av = np.asarray(Parallel(n_jobs=numb_jobs)(delayed(realization)(i, L) for i in range(numb_itr)))
    # for itr in range(numb_itr):
    #     ax0.plot(W_list, S_av[itr], '-', lw=0.1)
    ax0.plot(W_list, np.mean(S_av, axis=0), '.-', label=f'$L={L}$')

ax0.legend()
ax0.set_xlabel('$W$')
ax0.xaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax0.set_ylabel('$S$')
ax0.yaxis.set_major_formatter(FormatStrFormatter('$%g$'))
ax0.set_title(f'Ergodic to MBL phase transition ({numb_itr} disorders)')

plt.savefig(f"/home/bart/Documents/papers/MBF/Ponte_2015/entropy_scaling/entropy_W_scaling_Hamiltonian_L_{L}.png", bbox_inches='tight', dpi=300)
plt.show()
