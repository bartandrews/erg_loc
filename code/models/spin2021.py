# --- python imports
import numpy as np
# --- QuSpin imports
from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d


def _drive(t, Omega):
    return np.sign(np.cos(Omega*t))


def spin2021(L, Nup, pauli, J_1, J_2, T_0=10):

    Omega = 2 * np.pi / T_0
    drive_args = [Omega]

    basis = spin_basis_1d(L, a=2, Nup=Nup, pauli=pauli)

    G_1 = [[J_1, i] for i in range(L)]
    J_1_x = [[-J_1, i, (i + 1) % L] for i in range(0, L, 2)]
    J_1_y = [[-J_1, i, (i + 1) % L] for i in range(0, L, 2)]
    J_1_z = [[-J_1, i, (i + 1) % L] for i in range(0, L, 2)]
    pos_G_2 = [[J_2, i] for i in range(L)]
    pos_J_2_x = [[-J_2, i, (i - 1) % L] for i in range(0, L, 2)]
    pos_J_2_y = [[-J_2, i, (i - 1) % L] for i in range(0, L, 2)]
    pos_J_2_z = [[-J_2, i, (i - 1) % L] for i in range(0, L, 2)]
    neg_G_2 = [[-J_2, i] for i in range(L)]
    neg_J_2_x = [[J_2, i, (i - 1) % L] for i in range(0, L, 2)]
    neg_J_2_y = [[J_2, i, (i - 1) % L] for i in range(0, L, 2)]
    neg_J_2_z = [[J_2, i, (i - 1) % L] for i in range(0, L, 2)]
    static = [["I", G_1], ["xx", J_1_x], ["yy", J_1_y], ["zz", J_1_z],
              ["I", pos_G_2], ["xx", pos_J_2_x], ["yy", pos_J_2_y], ["zz", pos_J_2_z]]
    dynamic = [["I", G_1, _drive, drive_args], ["xx", J_1_x, _drive, drive_args],
               ["yy", J_1_y, _drive, drive_args], ["zz", J_1_z, _drive, drive_args],
               ["I", neg_G_2, _drive, drive_args], ["xx", neg_J_2_x, _drive, drive_args],
               ["yy", neg_J_2_y, _drive, drive_args], ["zz", neg_J_2_z, _drive, drive_args]]
    H = 0.5 * Omega * hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False,
                                  check_pcon=False)

    return H


def spin2021_2(L, Nup, pauli, J_1, J_2, W, h_0=2):

    basis = spin_basis_1d(L, Nup=Nup, pauli=pauli)

    # h_term = [[+h_0, L//2]]
    h_term = [[np.random.uniform(-W, W), i] for i in range(L)]
    static = [["z", h_term]]
    dynamic = []
    V = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False,
                    check_pcon=False)

    G_1 = [[J_1, i] for i in range(L)]
    J_1_x = [[-J_1, i, (i+1) % L] for i in range(0, L, 2)]
    J_1_y = [[-J_1, i, (i+1) % L] for i in range(0, L, 2)]
    J_1_z = [[-J_1, i, (i+1) % L] for i in range(0, L, 2)]
    # h_z_term = [[np.random.uniform(-W, W), i] for i in range(L)]
    static = [["I", G_1], ["xx", J_1_x], ["yy", J_1_y], ["zz", J_1_z]]
    dynamic = []
    H_1 = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False,
                      check_pcon=False)

    G_2 = [[J_2, i] for i in range(L)]
    J_2_x = [[-J_2, i, (i-1) % L] for i in range(0, L, 2)]
    J_2_y = [[-J_2, i, (i-1) % L] for i in range(0, L, 2)]
    J_2_z = [[-J_2, i, (i-1) % L] for i in range(0, L, 2)]
    # h_z_term = [[np.random.uniform(-W, W), i] for i in range(L)]
    static = [["I", G_2], ["xx", J_2_x], ["yy", J_2_y], ["zz", J_2_z]]
    dynamic = []
    H_2 = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False,
                      check_pcon=False)

    return V, H_1, H_2
