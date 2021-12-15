# --- python imports
import numpy as np
# --- QuSpin imports
from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d


def _V_drive(t, T_0, T_1):
    if 0 <= t < T_1 / 2:
        ans = 1
    elif T_1 / 2 <= t < T_0 + T_1 / 2:
        ans = 0
    elif T_1 / 2 + T_0 < T_0 + T_1:
        ans = 1
    else:
        ans = 0

    return (T_0 / np.pi) * ans


def _H1_drive(t, T_0, T_1):

    if 0 <= t < T_1/2:
        ans = 0
    elif T_1/2 <= t < T_1/2 + T_0/4:
        ans = 1
    elif T_1/2 + T_0/4 <= t < T_1/2 + 3*T_0/4:
        ans = 0
    elif T_1/2 + 3*T_0/4 <= t < T_0 + T_1/2:
        ans = 1
    else:
        ans = 0

    return ans


def _H2_drive(t, T_0, T_1):

    if 0 <= t < T_1/2:
        ans = 0
    elif T_1/2 <= t < T_1/2 + T_0/4:
        ans = 0
    elif T_1/2 + T_0/4 <= t < T_1/2 + 3*T_0/4:
        ans = 1
    elif T_1/2 + 3*T_0/4 <= t < T_0 + T_1/2:
        ans = 0
    else:
        ans = 0

    return ans


def spin2021(L, Nup, pauli, J_1, J_2, W, T_0=1, T_1=1):

    drive_args = [T_0, T_1]

    basis = spin_basis_1d(L, Nup=Nup, pauli=pauli)

    G_1 = [[J_1, i] for i in range(L)]
    J_1_x = [[-J_1, i, (i+1)] for i in range(0, L, 2)]
    J_1_y = [[-J_1, i, (i+1)] for i in range(0, L, 2)]
    J_1_z = [[-J_1, i, (i+1)] for i in range(0, L, 2)]
    G_2 = [[J_2, i] for i in range(L)]
    J_2_x = [[-J_2, i, (i-1)] for i in range(2, L, 2)]
    J_2_y = [[-J_2, i, (i-1)] for i in range(2, L, 2)]
    J_2_z = [[-J_2, i, (i-1)] for i in range(2, L, 2)]
    h_term = [[np.random.uniform(-W * np.pi, W * np.pi), i] for i in range(L)]
    static = []
    dynamic = [["I", G_1, _H1_drive, drive_args], ["xx", J_1_x, _H1_drive, drive_args], ["yy", J_1_y, _H1_drive, drive_args], ["zz", J_1_z, _H1_drive, drive_args],
               ["I", G_2, _H2_drive, drive_args], ["xx", J_2_x, _H2_drive, drive_args], ["yy", J_2_y, _H2_drive, drive_args], ["zz", J_2_z, _H2_drive, drive_args],
               ["z", h_term, _V_drive, drive_args]]
    H = (np.pi/T_0) * hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False,
                                  check_pcon=False)

    return H


def spin2021_2(L, Nup, pauli, J_1, J_2, W):

    basis = spin_basis_1d(L, Nup=Nup, pauli=pauli)

    h_term = [[np.random.uniform(-W*np.pi, W*np.pi), i] for i in range(L)]
    static = [["z", h_term]]
    dynamic = []
    V = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False,
                    check_pcon=False)

    G_1 = [[J_1, i] for i in range(L)]
    J_1_x = [[-J_1, i, (i+1)] for i in range(0, L, 2)]
    J_1_y = [[-J_1, i, (i+1)] for i in range(0, L, 2)]
    J_1_z = [[-J_1, i, (i+1)] for i in range(0, L, 2)]
    static = [["I", G_1], ["xx", J_1_x], ["yy", J_1_y], ["zz", J_1_z]]
    dynamic = []
    H_1 = np.pi * hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False,
                              check_pcon=False)

    G_2 = [[J_2, i] for i in range(L)]
    J_2_x = [[-J_2, i, (i-1)] for i in range(2, L, 2)]
    J_2_y = [[-J_2, i, (i-1)] for i in range(2, L, 2)]
    J_2_z = [[-J_2, i, (i-1)] for i in range(2, L, 2)]
    static = [["I", G_2], ["xx", J_2_x], ["yy", J_2_y], ["zz", J_2_z]]
    dynamic = []
    H_2 = np.pi * hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False,
                              check_pcon=False)

    return V, H_1, H_2
