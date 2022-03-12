# --- python imports
import numpy as np
# --- QuSpin imports
from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian


def pxp(L, Nup, pauli, bc, J):

    # pblock=1, a=1, kblock=0
    basis = spin_basis_1d(L, Nup=Nup, pauli=pauli)

    if bc is "p":
        first_term = [[J/4, i, (i+1) % L, (i+2) % L] for i in range(L)]
        second_term = [[-J/4, i, (i+1) % L, (i+2) % L] for i in range(L)]
        third_term = [[-J/4, i, (i+1) % L, (i+2) % L] for i in range(L)]
        fourth_term = [[J/4, i, (i+1) % L, (i+2) % L] for i in range(L)]
        static = [["IxI", first_term], ["Ixz", second_term], ["zxI", third_term], ["zxz", fourth_term]]
        dynamic = []
    else:  # "o"
        first_term = [[J/4, i] for i in range(L)]
        second_term = [[-J/4, i, i+1] for i in range(L-1)]
        third_term = [[-J/4, i-1, i] for i in range(1, L)]
        fourth_term = [[J/4, i-1, i, i+1] for i in range(1, L-1)]
        initial_term = [[J/2, 0, 1]]
        minus_initial_term = [[-J/2, 0, 1]]
        final_term = [[J/2, L-2, L-1]]
        minus_final_term = [[-J/2, L-2, L-1]]
        static = [["x", first_term], ["xz", second_term], ["zx", third_term], ["zxz", fourth_term],
                  ["xI", initial_term], ["xz", minus_initial_term], ["Ix", final_term], ["zx", minus_final_term]]
        dynamic = []

    H = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False,
                    check_pcon=False)

    return H
