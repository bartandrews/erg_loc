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

    leaf = fp.file_name_leaf("W_flow", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("W_flow", path, model, leaf)

    # "ener_W_flow", "r_W_flow", "ent_W_flow"
    tools = ["ener_W_flow", "r_W_flow", "ent_W_flow"]
    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    def realization(itr, _model, _leaf_args, entropy=False):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        _E_array = np.zeros(_leaf_args['W_samp'], dtype=object)
        _S_array = np.zeros(_leaf_args['W_samp'], dtype=object)

        _r = []
        for i, _W in enumerate(np.linspace(_leaf_args['W_min'], _leaf_args['W_max'], _leaf_args['W_samp'])):
            _leaf_args['W'] = _W
            H = fh.chosen_hamiltonian(_model, _leaf_args)
            if entropy:
                _E, psi = H.eigh()
                _E_array[i] = _E

                r_tmp = []
                for j in range(1, len(_E)-1):
                    delta_n = _E[j] - _E[j-1]
                    delta_n_plus_1 = _E[j+1] - _E[j]
                    r_tmp.append(min(delta_n, delta_n_plus_1) / max(delta_n, delta_n_plus_1))
                _r.append(np.mean(r_tmp))

                _S_array[i] = []
                for j in range(psi.shape[1]):
                    _S_array[i].append(float(H.basis.ent_entropy(psi[:, j], sub_sys_A=range(H.basis.L//2))["Sent_A"]))
                _S_array[i] = np.asarray(_S_array[i])
            else:
                _E = H.eigvalsh()
                _E_array[i] = _E

                r_tmp = []
                for j in range(1, len(_E) - 1):
                    delta_n = _E[j] - _E[j - 1]
                    delta_n_plus_1 = _E[j + 1] - _E[j]
                    r_tmp.append(min(delta_n, delta_n_plus_1) / max(delta_n, delta_n_plus_1))
                _r.append(np.mean(r_tmp))

        if entropy:
            return _E_array, _r, _S_array
        else:
            return _E_array, _r, None

    ###################################################################################################################

    ent_flag = True if "ent_W_flow" in tools else False
    array = np.asarray(Parallel(n_jobs=threads)(delayed(realization)(i, model, leaf_args, entropy=ent_flag)
                                                for i in range(leaf_args['dis'])), dtype=object)

    W_list = np.linspace(_leaf_args['W_min'], _leaf_args['W_max'], _leaf_args['W_samp'])
    W_ext_list = []
    for W in W_list:
        W_ext_list += [W]*len(array[0, 0, 0])

    ###############
    # ener_W_flow #
    ###############

    if "ener_W_flow" in tools:

        E_array = array[:, 0, :]
        E = np.mean(E_array, axis=0)
        E_list = [i.tolist() for i in E]
        E_list_flat = [item for sublist in E_list for item in sublist]

        for i, E_val in enumerate(E_list_flat):
            data['ener_W_flow'].write(f"{W_ext_list[i]:g}\t{E_val}\n")

    ############
    # r_W_flow #
    ############

    if "r_W_flow" in tools:

        r_array = array[:, 1, :]
        r = np.mean(r_array, axis=0)

        for i, r_val in enumerate(r):
            data['r_W_flow'].write(f"{W_list[i]:g}\t{r_val}\n")

    ##############
    # ent_W_flow #
    ##############

    if "ent_W_flow" in tools:

        S_array = array[:, 2, :]
        S = np.mean(S_array, axis=0)
        S_list = [i.tolist() for i in S]
        S_list_flat = [item for sublist in S_list for item in sublist]

        for i, S_val in enumerate(S_list_flat):
            data['ent_W_flow'].write(f"{W_ext_list[i]:g}\t{S_val}\n")

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('W_flow')

    my_W_flow(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
