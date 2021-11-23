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
W = 0.5  # W=0.5 for ergodic and W=8 for MBL
# V parameters
h_0 = 2
# driving parameters
T_0 = 7
T_1_max = 3
T_1_samp = 16
# iteration parameters
numb_itr = 10  # 20000 for L=8,10 or 1000 for L=12,14
# initialize results arrays
PR_sum = np.zeros((numb_itr, T_1_samp))
PR = np.zeros((numb_itr, T_1_samp))

t_0 = time.time()
#
for itr in range(numb_itr):
    print(f"Iteration {itr+1} of {numb_itr}")
    for i, T_1 in enumerate(np.linspace(0, T_1_max, T_1_samp)):
        # print(f"T_1 = {T_1}")
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

        for alpha_idx in range(len(alpha)):
            PR_sum[itr][i] += np.abs(np.dot(psi[:, 0], alpha[:, alpha_idx]))**4

        PR[itr][i] = 1/(2**L * PR_sum[itr][i])
#
print(f"Total time (seconds): {int(time.time()-t_0)}")
print(f"Average time per iteration (seconds): {(time.time()-t_0)/numb_itr:.1f}")


fig = plt.figure()
ax = fig.add_subplot(111)
T_1s = np.linspace(0, T_1_max, T_1_samp)
# plot the iterations
for itr in range(numb_itr):
    ax.plot(T_1s, PR[itr], '-', lw=0.1)
# plot the average
ax.plot(T_1s, np.mean(PR, axis=0), '-', marker='x', c='k', lw=1)
ax.set_xlabel("$T_1$")
ax.set_ylabel("PR")
ax.set_title(f"$W={W}$, $L={L}$")
plt.show()
