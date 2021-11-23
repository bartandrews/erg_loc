# --- python imports
import numpy as np
# --- QuSpin imports
from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d


def heisenberg(L, Nup, pauli, J_x, J_y, J_z, W):

    basis = spin_basis_1d(L, Nup=Nup, pauli=pauli)

    J_x_term = [[J_x, i, i+1] for i in range(L-1)]
    J_y_term = [[J_y, i, i+1] for i in range(L-1)]
    J_z_term = [[J_z, i, i+1] for i in range(L-1)]
    h_z_term = [[np.random.uniform(-W, W), i] for i in range(L)]
    static = [["xx", J_x_term], ["yy", J_y_term], ["zz", J_z_term], ["z", h_z_term]]
    dynamic = []

    H = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False,
                    check_pcon=False)

    return H
