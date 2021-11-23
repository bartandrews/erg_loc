import numpy as np
import matplotlib.pyplot as plt
import time

from quspin.operators import hamiltonian, exp_op
from quspin.basis import spin_basis_1d

# system parameters
L = 8
basis = spin_basis_1d(L)
# H_0 parameters
J_x_0 = 1
J_z_0 = 1
W = 8  # W=0.5 for ergodic and W=8 for MBL
# V parameters
h_0 = 2
# driving parameters
T_0 = 7
T_1 = 1.5


t_0 = time.time()
#
# compute H_0
J_x = [[J_x_0, i, i+1] for i in range(L-1)]
J_y = [[J_x_0, i, i+1] for i in range(L-1)]
J_z = [[J_z_0, i, i+1] for i in range(L-1)]
h_z = [[np.random.uniform(-W, W), i] for i in range(L)]
static_H_0 = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z]]
dynamic_H_0 = []
H_0 = hamiltonian(static_H_0, dynamic_H_0, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)

# compute V
h = [[h_0, L//2]]
static_V = [["z", h]]
dynamic_V = []
V = hamiltonian(static_V, dynamic_V, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)

# --- eigenstate overlap magnitude
E_Tinf = H_0.trace() / (2**L)
E, alpha = H_0.eigh()

F = exp_op(H_0*T_0+V*T_1, a=-1j)
# F_H = exp_op(H_0*T_0, a=-1j)
# F_V = exp_op(V*T_1, a=-1j)
# F = F_H.dot(F_V.dot(np.ones(len(F_V.get_mat(dense=True)))))

epsilon, psi = np.linalg.eig(F.get_mat(dense=True))

A2 = np.zeros((len(psi), len(alpha)))
for i_idx in range(len(psi)):
    for alpha_idx in range(len(alpha)):
        A2[alpha_idx, i_idx] = np.abs(np.dot(psi[:, i_idx], alpha[:, alpha_idx]))**2
#
print(f"Total time (seconds): {time.time()-t_0:.1f}")


fig = plt.figure()
ax = fig.add_subplot(111)
for i_idx in range(len(psi)):
    ax.plot(range(len(alpha)), A2[:, i_idx], '.', label=f"$i={i_idx}$")
ax.set_xlabel("$\\alpha$")
ax.set_ylabel("$|A_{\\alpha, i}|^2$")
ax.set_yscale('log')
ax.set_ylim([10e-12, 0.01])
ax.set_title(f"$W={W}$, $L={L}$")
ax.axhline(0.0001, c='r')
plt.show()
