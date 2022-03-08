# --- python imports
import numpy as np
from time import perf_counter
import sys
from joblib import delayed, Parallel
# --- driven_systems imports
import functions.func_ham as fh
import functions.func_args as fa
import functions.func_proc as fp


def find_Ns(_model, _leaf_args):

    _leaf_args['W'] = _leaf_args['W_min']
    H = fh.chosen_hamiltonian(_model, _leaf_args)

    return H.Ns


def my_W_flow(path_flag, threads, model, _leaf_args):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = perf_counter()  # start the timer

    leaf = fp.file_name_leaf("W_flow", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("W_flow", path, model, leaf)

    # "ener_W_flow", "r_W_flow", "ent_W_flow", "ent_rel_W_flow"
    tools = ["ent_rel_W_flow"]
    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    def realization(itr, _model, _leaf_args, entropy=False):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        Ns = find_Ns(_model, _leaf_args)

        _ener_array = np.zeros((_leaf_args['W_samp'], Ns))
        _r_array = np.zeros((_leaf_args['W_samp'], Ns))
        _ent_array = np.zeros((_leaf_args['W_samp'], Ns))
        _ent_rel_array = np.zeros((_leaf_args['W_samp'], Ns))

        for i, _W in enumerate(np.linspace(_leaf_args['W_min'], _leaf_args['W_max'], _leaf_args['W_samp'])):
            _leaf_args['W'] = _W
            H = fh.chosen_hamiltonian(_model, _leaf_args)
            if entropy:
                E, psi = H.eigh()

                if "ener_W_flow" in tools:
                    for j, ener_val in enumerate(E):
                        _ener_array[i, j] = ener_val

                if "r_W_flow" in tools:
                    _r_array[i, 0] = None
                    _r_array[i, Ns-1] = None
                    for j in range(1, Ns-1):
                        delta_n = E[j] - E[j-1]
                        delta_n_plus_1 = E[j+1] - E[j]
                        _r_array[i, j] = min(delta_n, delta_n_plus_1) / max(delta_n, delta_n_plus_1)

                if "ent_W_flow" in tools:
                    for j in range(Ns):
                        _ent_array[i, j] = float(H.basis.ent_entropy(psi[:, j],
                                                                     sub_sys_A=range(H.basis.L//2))["Sent_A"])

                if "ent_rel_W_flow" in tools:
                    E_min, E_max = H.eigsh(k=2, which="BE", maxiter=1E4, return_eigenvectors=False)
                    E_target = E_min + 0.5 * (E_max - E_min)
                    _, psi = H.eigsh(k=2, sigma=E_target, maxiter=1E4)
                    p = np.empty(Ns, dtype=complex)
                    q = np.empty(Ns, dtype=complex)
                    H.basis.inplace_Op(psi[:, 0], [["z", [i], 1] for i in range(H.basis.L)], dtype=complex, v_out=p)
                    H.basis.inplace_Op(psi[:, 1], [["z", [i], 1] for i in range(H.basis.L)], dtype=complex, v_out=q)
                    for j in range(Ns):
                        _ent_rel_array[i, j] = np.abs(psi[j, 0])**2 * 2 * np.log(np.abs(psi[j, 0]/psi[j, 1]))
            else:
                E = H.eigvalsh()

                if "ener_W_flow" in tools:
                    for j, ener_val in enumerate(E):
                        _ener_array[i, j] = ener_val

                if "r_W_flow" in tools:
                    _r_array[i, 0] = None
                    _r_array[i, Ns-1] = None
                    for j in range(1, Ns-1):
                        delta_n = E[j] - E[j-1]
                        delta_n_plus_1 = E[j+1] - E[j]
                        _r_array[i, j] = min(delta_n, delta_n_plus_1) / max(delta_n, delta_n_plus_1)

        if entropy:
            return _ener_array, _r_array, _ent_array, _ent_rel_array
        else:
            return _ener_array, _r_array

    ###################################################################################################################

    if any(item in ["r_W_flow", "ent_W_flow", "ent_rel_W_flow"] for item in tools):
        ent_flag = True
    else:
        ent_flag = False

    W_list = np.linspace(_leaf_args['W_min'], _leaf_args['W_max'], _leaf_args['W_samp'])
    array = np.stack(Parallel(n_jobs=threads)(delayed(realization)(i, model, leaf_args, entropy=ent_flag)
                                              for i in range(leaf_args['dis'])), axis=0)
    # (disorder, tool, samp, state)

    ###############
    # ener_W_flow #
    ###############

    if "ener_W_flow" in tools:

        ener_array = array[:, 0]
        ener = np.mean(ener_array, axis=0)

        for i, W_val in enumerate(W_list):
            for ener_state_val in ener[i, :]:
                data['ener_W_flow'].write(f"{W_val:g}\t{ener_state_val}\n")

    ############
    # r_W_flow #
    ############

    if "r_W_flow" in tools:

        r_array = array[:, 1, :, 1:-1]  # truncate first and last None values
        r = np.mean(r_array, axis=(0, 2))

        for i, r_val in enumerate(r):
            data['r_W_flow'].write(f"{W_list[i]:g}\t{r_val}\n")

    ##############
    # ent_W_flow #
    ##############

    if "ent_W_flow" in tools:

        ent_array = array[:, 2]
        ent = np.mean(ent_array, axis=0)

        for i, W_val in enumerate(W_list):
            for ent_state_val in ent[i, :]:
                data['ent_W_flow'].write(f"{W_val:g}\t{ent_state_val}\n")

    ##################
    # ent_rel_W_flow #
    ##################

    if "ent_rel_W_flow" in tools:

        ent_rel_array = np.sum(array[:, 3], axis=2)
        ent_rel = np.mean(ent_rel_array, axis=0)

        for i, ent_rel_val in enumerate(ent_rel):
            data['ent_rel_W_flow'].write(f"{W_list[i]:g}\t{ent_rel_val}\n")

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('W_flow')

    my_W_flow(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
