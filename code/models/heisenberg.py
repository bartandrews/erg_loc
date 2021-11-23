import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
from matplotlib import ticker
import time
from joblib import delayed, Parallel

from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d


class HeisenbergModel(hamiltonian):

    def __init__(self, params):
        hamiltonian.__init__(self, params)

    def init_model(self, params):

        basis = spin_basis_1d(L)

        J_x = [[J_x_0, i, i + 1] for i in range(L - 1)]
        J_y = [[J_x_0, i, i + 1] for i in range(L - 1)]
        J_z = [[J_z_0, i, i + 1] for i in range(L - 1)]
        h_z = [[np.random.uniform(-W_val, W_val), i] for i in range(L)]
        static = [["xx", J_x], ["yy", J_y], ["zz", J_z], ["z", h_z]]
        dynamic = []
        H = hamiltonian(static, dynamic, basis=basis, dtype=np.float64, check_symm=False, check_herm=False)
