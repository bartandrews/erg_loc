import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib import ticker

from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d
from quspin.tools.Floquet import Floquet

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


# compute H_1
J_x = [[J_x_0, i, i+1] for i in range(L-1)]
J_y = [[J_x_0, i, i+1] for i in range(L-1)]
J_z = [[J_z_0, i, i+1] for i in range(L-1)]
h_z = [[np.random.uniform(-0.5, 0.5), i] for i in range(L)]
pos_h = [[+h_0, L//2]]
neg_h = [[-h_0, L//2]]
static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z], ["z", pos_h]]
dynamic = [["xx", J_x, drive, drive_args], ["yy", J_y, drive, drive_args], ["zz", J_z, drive, drive_args],
           ["z", h_z, drive, drive_args], ["z", neg_h, drive, drive_args]]
H_1 = 0.5*hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)
# eigenbasis of H_1
E_1, alpha_1 = H_1.eigh(time=0)
# eigenbasis of Floq_1
t_list = np.array([0.0, T_0/2.0, T_0/2.0+T_1]) + np.finfo(float).eps
dt_list = np.array([T_0/2.0, T_1, T_0/2.0])
Floq_1 = Floquet({'H': H_1, 't_list': t_list, 'dt_list': dt_list}, VF=True)
psi_1 = Floq_1.VF
# compute A2_1
A2_1 = np.zeros((len(psi_1), len(alpha_1)))
for i_idx in range(len(psi_1)):
    for alpha_idx in range(len(alpha_1)):
        A2_1[alpha_idx, i_idx] = np.abs(np.dot(psi_1[:, i_idx], alpha_1[:, alpha_idx]))**2

# compute H_2
J_x = [[J_x_0, i, i + 1] for i in range(L - 1)]
J_y = [[J_x_0, i, i + 1] for i in range(L - 1)]
J_z = [[J_z_0, i, i + 1] for i in range(L - 1)]
h_z = [[np.random.uniform(-8, 8), i] for i in range(L)]
pos_h = [[+h_0, L // 2]]
neg_h = [[-h_0, L // 2]]
static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z], ["z", pos_h]]
dynamic = [["xx", J_x, drive, drive_args], ["yy", J_y, drive, drive_args], ["zz", J_z, drive, drive_args],
           ["z", h_z, drive, drive_args], ["z", neg_h, drive, drive_args]]
H_2 = 0.5 * hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)
# eigenbasis of H_2
E_2, alpha_2 = H_2.eigh(time=0)
# eigenbasis of Floq_2
t_list = np.array([0.0, T_0 / 2.0, T_0 / 2.0 + T_1]) + np.finfo(float).eps
dt_list = np.array([T_0 / 2.0, T_1, T_0 / 2.0])
Floq_2 = Floquet({'H': H_2, 't_list': t_list, 'dt_list': dt_list}, VF=True)
psi_2 = Floq_2.VF
# compute A2_2
A2_2 = np.zeros((len(psi_2), len(alpha_2)))
for i_idx in range(len(psi_2)):
    for alpha_idx in range(len(alpha_2)):
        A2_2[alpha_idx, i_idx] = np.abs(np.dot(psi_2[:, i_idx], alpha_2[:, alpha_idx])) ** 2

fig = plt.figure(figsize=(10, 5))
fig.suptitle("Floquet Eigenstate Structure, $\ket{\phi_N}=\sum_i A_{\\alpha_0, i} e^{-\mathrm{i}\omega_i NT}\ket{\psi_i} \;\Rightarrow\; A_{\\alpha,i}=\\braket{\psi_i|\\alpha}$")
gs = gridspec.GridSpec(1, 2, hspace=0, wspace=0)
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1], sharey=ax0)
# plot ergodic
ax0.plot(range(len(alpha_1)), A2_1[:, 0], '.')
ax0.set_xlabel("$\\alpha$")
ax0.set_ylabel("$|A_{\\alpha, 0}|^2$")
ax0.set_yscale('log')
ax0.set_ylim([10e-15, 1])
ax0.set_title(f"Ergodic ($W=0.5$, $L={L}$)")
ax0.axhline(0.0001, c='r')
# plot MBL
ax1.plot(range(len(alpha_2)), A2_2[:, 0], '.')
ax1.set_xlabel("$\\alpha$")
ax1.yaxis.set_visible(False)
ax1.set_title(f"MBL ($W=8$, $L={L}$)")
ax1.axhline(0.0001, c='r')

plt.savefig("/home/bart/Documents/papers/MBF/Ponte_2015/floquet_eigenstate_structure.png", bbox_inches='tight', dpi=300)
plt.show()
