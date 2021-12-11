# --- driven_systems imports
from models.heisenberg import heisenberg
from models.ponte2015 import ponte2015, ponte2015_2
from models.spin2021 import spin2021, spin2021_2


def chosen_hamiltonian(model, leaf_args):

    if model == "heisenberg":
        H = heisenberg(leaf_args['L'], leaf_args['Nup'], leaf_args['pauli'],
                       leaf_args['J'][0], leaf_args['J'][1], leaf_args['J'][2], leaf_args['W'])
    elif model == "ponte2015":
        H = ponte2015(leaf_args['L'], leaf_args['Nup'], leaf_args['pauli'],
                      leaf_args['J'][0], leaf_args['J'][1], leaf_args['J'][2], leaf_args['W'],
                      leaf_args['h0'], leaf_args['T0'], leaf_args['T1'])
    elif model == "ponte2015_2":
        H_0, V = ponte2015_2(leaf_args['L'], leaf_args['Nup'], leaf_args['pauli'],
                             leaf_args['J'][0], leaf_args['J'][1], leaf_args['J'][2], leaf_args['W'],
                             leaf_args['h0'])
    elif model == "spin2021":
        H = spin2021(leaf_args['L'], leaf_args['Nup'], leaf_args['pauli'],
                               leaf_args['J'][0], leaf_args['J'][1], leaf_args['W'])
    elif model == "spin2021_2":
        V, H_1, H_2 = spin2021_2(leaf_args['L'], leaf_args['Nup'], leaf_args['pauli'],
                               leaf_args['J'][0], leaf_args['J'][1], leaf_args['W'],
                               leaf_args['h0'])
    else:
        raise ValueError("model not implemented in func_ham")

    if model == "ponte2015_2":
        return H_0, V
    if model == "spin2021_2":
        return V, H_1, H_2
    else:
        return H
