# --- python imports
import numpy as np
# --- QuSpin imports
from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d


def drive(t, h_0, W, omega):
    val = h_0 + W*np.sin(omega*t)
    return val


def ising_2(L, Nup, pauli, J_x, J_y, J_z, W=1, h_0=2.3, T_0=1.5708):

    omega = 2*np.pi/T_0
    drive_args = [h_0, W, omega]

    basis = spin_basis_1d(L, Nup=Nup, pauli=pauli)

    J_z_term = [[J_z, i, i+1] for i in range(L-1)]
    J_y_term = [[J_y, i, i+1] for i in range(L-1)]
    J_x_term = [[J_x, i] for i in range(L)]

    static = [["zz", J_z_term], ["yy", J_y_term]]
    dynamic = [["x", J_x_term, drive, drive_args]]

    H = -0.5 * hamiltonian(static, dynamic, basis=basis, dtype=np.float64,
                           check_symm=False, check_herm=False, check_pcon=False)

    return H
