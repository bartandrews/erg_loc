# --- driven_systems imports
from models.heisenberg import heisenberg
from models.ponte2015 import ponte2015


def chosen_hamiltonian(model, leaf_args):

    if model == "heisenberg":
        H = heisenberg(leaf_args['L'], leaf_args['Nup'], leaf_args['pauli'],
                       leaf_args['J'][0], leaf_args['J'][1], leaf_args['J'][2], leaf_args['W'])
    elif model == "ponte2015":
        H = ponte2015(leaf_args['L'], leaf_args['Nup'], leaf_args['pauli'],
                      leaf_args['J'][0], leaf_args['J'][1], leaf_args['J'][2], leaf_args['W'],
                      leaf_args['h0'], leaf_args['T0'], leaf_args['T1'])
    else:
        raise ValueError("model not implemented in func_ham")

    return H
