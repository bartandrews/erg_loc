# --- python imports
import numpy as np
from time import perf_counter
import sys
from joblib import delayed, Parallel
# --- driven_systems imports
import functions.func_ham as fh
import functions.func_args as fa
import functions.func_proc as fp


def my_W_flow(path_flag, threads, model, _leaf_args):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = perf_counter()  # start the timer

    leaf = fp.file_name_leaf("L_flow", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("L_flow", path, model, leaf)

    tools = ["ener_L_flow", "ent_L_flow"]
    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    def realization(itr, _model, _leaf_args, entropy=False):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        _E_array = np.zeros(_leaf_args['L_samp'], dtype=object)
        _S_array = np.zeros(_leaf_args['L_samp'], dtype=object)

        _L_list = list(map(int, np.linspace(_leaf_args['L_min'], _leaf_args['L_max'], _leaf_args['L_samp'])))
        for i, _L in enumerate(_L_list):
            _leaf_args['L'] = _L
            H = fh.chosen_hamiltonian(_model, _leaf_args)
            if entropy:
                _E, psi = H.eigh()
                _E_array[i] = _E
                _S_array[i] = []
                for j in range(psi.shape[1]):
                    _S_array[i].append(float(H.basis.ent_entropy(psi[:, j], sub_sys_A=range(H.basis.L//2))["Sent_A"]))
                _S_array[i] = np.asarray(_S_array[i])
            else:
                _E = H.eigvalsh()
                _E_array[i] = _E

        if entropy:
            return _E_array, _S_array
        else:
            return _E_array, None

    ###################################################################################################################

    ent_flag = True if "ent_L_flow" in tools else False
    array = np.asarray(Parallel(n_jobs=threads)(delayed(realization)(i, model, leaf_args, entropy=ent_flag)
                                                for i in range(leaf_args['dis'])), dtype=object)

    L_list = list(map(int, np.linspace(_leaf_args['L_min'], _leaf_args['L_max'], _leaf_args['L_samp'])))
    L_ext_list = []
    for i, L in enumerate(L_list):
        L_ext_list += [L]*len(array[0, 0, i])

    ###############
    # ener_L_flow #
    ###############

    if "ener_L_flow" in tools:

        E_array = array[:, 0, :]
        E = np.mean(E_array, axis=0)
        E_list = [i.tolist() for i in E]
        E_list_flat = [item for sublist in E_list for item in sublist]

        for i, E_val in enumerate(E_list_flat):
            data['ener_L_flow'].write(f"{L_ext_list[i]}\t{E_val}\n")

    ##############
    # ent_L_flow #
    ##############

    if "ent_L_flow" in tools:

        S_array = array[:, 1, :]
        S = np.mean(S_array, axis=0)
        S_list = [i.tolist() for i in S]
        S_list_flat = [item for sublist in S_list for item in sublist]

        for i, S_val in enumerate(S_list_flat):
            data['ent_L_flow'].write(f"{L_ext_list[i]}\t{S_val}\n")

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('L_flow')

    my_W_flow(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
