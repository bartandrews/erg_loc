# --- driven_systems imports
from models.heisenberg import heisenberg


def chosen_hamiltonian(model, leaf_args):

    if model == "heisenberg":
        H = heisenberg(leaf_args['L'], leaf_args['Nup'], leaf_args['pauli'],
                       leaf_args['J'][0], leaf_args['J'][1], leaf_args['J'][2], leaf_args['W'])
    else:
        raise ValueError("model not implemented in func_ham")

    return H
