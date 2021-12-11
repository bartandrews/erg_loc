# --- python imports
import numpy as np
# --- QuSpin imports
from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d


def _drive(t, T_0, T_1):

    T = T_0 + T_1

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


def ponte2015(L, Nup, pauli, J_x, J_y, J_z, W, h_0=2, T_0=7, T_1=1.5):

    basis = spin_basis_1d(L, Nup=Nup, pauli=pauli)

    J_x_term = [[J_x, i, i+1] for i in range(L-1)]
    J_y_term = [[J_y, i, i+1] for i in range(L-1)]
    J_z_term = [[J_z, i, i+1] for i in range(L-1)]
    h_z_term = [[np.random.uniform(-W, W), i] for i in range(L)]
    pos_h_term = [[+h_0, L//2]]
    neg_h_term = [[-h_0, L//2]]
    static = [["xx", J_x_term], ["yy", J_y_term], ["zz", J_z_term], ["z", h_z_term], ["z", pos_h_term]]
    dynamic = [["xx", J_x_term, _drive, [T_0, T_1]], ["yy", J_y_term, _drive, [T_0, T_1]],
               ["zz", J_z_term, _drive, [T_0, T_1]], ["z", h_z_term, _drive, [T_0, T_1]],
               ["z", neg_h_term, _drive, [T_0, T_1]]]
    H = 0.5 * hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False,
                          check_pcon=False)

    return H


def ponte2015_2(L, Nup, pauli, J_x, J_y, J_z, W, h_0=2):

    basis = spin_basis_1d(L, Nup=Nup, pauli=pauli)

    J_x_term = [[J_x, i, i + 1] for i in range(L - 1)]
    J_y_term = [[J_y, i, i + 1] for i in range(L - 1)]
    J_z_term = [[J_z, i, i + 1] for i in range(L - 1)]
    h_z_term = [[np.random.uniform(-W, W), i] for i in range(L)]
    static = [["xx", J_x_term], ["yy", J_y_term], ["zz", J_z_term], ["z", h_z_term]]
    dynamic = []
    H_0 = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False,
                      check_pcon=False)

    h_term = [[h_0, L//2]]
    static = [["z", h_term]]
    dynamic = []
    V = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False,
                    check_pcon=False)

    return H_0, V
