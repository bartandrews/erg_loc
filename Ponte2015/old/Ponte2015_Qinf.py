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
numb_N = 31
# initialize results arrays
Q_N = np.zeros((numb_itr, T_1_samp))

t_0 = time.time()
#
for itr in range(numb_itr):
    print(f"Iteration {itr + 1} of {numb_itr}")
    for i, T_1 in enumerate(np.linspace(0, T_1_max, T_1_samp)):

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

        # --- energy absorbed under driving
        E_Tinf = H_0.trace() / (2**L)
        E, phi = H_0.eigh()
        # print(f"E = {E}")
        # print(f"phi[:,0] = {phi[:,0]}")
        # print(f"phi[:,1] = {phi[:,1]}")
        # print(f"phi[:,2] = {phi[:,2]}")
        # print(f"phi[:,3] = {phi[:,3]}")
        # print(f"H_0.phi[:,0] = {H_0.dot(phi[:, 0])}")
        # print(f"H_0.phi[:,1] = {H_0.dot(phi[:, 1])}")
        # print(f"H_0.phi[:,2] = {H_0.dot(phi[:, 2])}")
        # print(f"H_0.phi[:,3] = {H_0.dot(phi[:, 3])}")
        E_0 = np.min(E)
        phi_0 = phi[:, np.argmin(E)]
        # print(f"E_0 = {E_0}")
        # print(f"phi_0 = {phi_0}")

        # F = exp_op(H_0*T_0+V*T_1, a=-1j)
        F_H = exp_op(H_0*T_0, a=-1j)
        F_V = exp_op(V*T_1, a=-1j)

        phi_N = phi_0
        for n in range(numb_N):
            # F_N = exp_op(H_0*T_0+V*T_1, a=-1j*n)
            # phi_N = F_N.dot(phi_0)

            if n > 0:
                phi_N = F_H.dot(F_V.dot(phi_N))

        Q_N[itr][i] = (np.real(H_0.matrix_ele(phi_N.conj().T, phi_N)) - E_0) / (E_Tinf - E_0)
#
print(f"Total time (seconds): {int(time.time()-t_0)}")
print(f"Average time per iteration (seconds): {(time.time()-t_0)/numb_itr:.1f}")


fig = plt.figure()
ax = fig.add_subplot(111)
T_1s = np.linspace(0, T_1_max, T_1_samp)
# plot the iterations
for itr in range(numb_itr):
    ax.plot(T_1s, Q_N[itr], '-', lw=0.1)
# plot the average
ax.plot(T_1s, np.mean(Q_N, axis=0), '-', marker='x', c='k', lw=1)
ax.set_xlabel("$T_1$")
ax.set_ylabel("$Q_\infty$")
ax.set_title(f"$W={W}$, $L={L}$, averaged over {numb_itr} disorders")
ax.axhline(1, c='r')
plt.show()
