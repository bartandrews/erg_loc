# --- python imports
import numpy as np
from time import perf_counter
import sys
from joblib import delayed, Parallel
# --- driven_systems imports
import functions.func_ham as fh
import functions.func_args as fa
import functions.func_proc as fp


def my_inst_ham(path_flag, threads, model, _leaf_args):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = perf_counter()  # start the timer

    leaf = fp.file_name_leaf("inst_ham", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("inst_ham", path, model, leaf)

    # "ener", "ener_spac", "ent", "ent_mid"
    tools = ["ener", "ener_spac", "ent"]
    if "ent_mid" in tools and len(tools) != 1:
        raise ValueError("The tool ent_mid can only be used in isolation.")

    if "ent_mid" not in tools:
        data = fp.prepare_output_files(tools, path, model, leaf)
    else:
        data = fp.prepare_output_files(["ent_mid"], path, model, leaf)

    ###################################################################################################################

    def realization(itr, _model, _leaf_args, entropy=0):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        H = fh.chosen_hamiltonian(_model, _leaf_args)

        _ener_array = np.zeros(H.Ns)
        _ener_spac_array = np.zeros(H.Ns)
        _ent_array = np.zeros(H.Ns)

        if entropy == 2:  # ent
            E, psi = H.eigh()

            # --- ener
            for i, ener_val in enumerate(E):
                _ener_array[i] = ener_val
            # --- ener_spac
            _ener_spac_array[H.Ns-1] = None
            for i in range(H.Ns-1):
                _ener_spac_array[i] = E[i+1]-E[i]
            # --- ent
            for i in range(H.Ns):
                _ent_array[i] = H.basis.ent_entropy(psi[:, i], sub_sys_A=range(H.basis.L//2))["Sent_A"]

            return _ener_array, _ener_spac_array, _ent_array
        elif entropy == 1:  # ent_mid
            _, psi = H.eigsh(k=1, sigma=0.5, maxiter=1E4)

            # --- ent_mid
            _ent_mid = H.basis.ent_entropy(psi, sub_sys_A=range(H.basis.L//2))["Sent_A"]

            return None, None, _ent_mid
        elif entropy == 0:  # no ent
            E = H.eigvalsh()

            # --- ener
            for i, ener_val in enumerate(E):
                _ener_array[i] = ener_val
            # --- ener_spac
            _ener_spac_array[H.Ns - 1] = None
            for i in range(H.Ns - 1):
                _ener_spac_array = E[i + 1] - E[i]

            return _ener_array, _ener_spac_array, None
        else:
            raise ValueError("ent_flag not defined in realization")

    ###################################################################################################################

    ent_flag = 0
    if any("ent" in i for i in tools):
        if "ent_mid" in tools:
            ent_flag = 1
        elif "ent" in tools:
            ent_flag = 2

    array = np.stack(Parallel(n_jobs=threads)(delayed(realization)(i, model, leaf_args, entropy=ent_flag)
                                              for i in range(leaf_args['dis'])), axis=0)  # (disorder, tool, state)

    ########
    # ener #
    ########

    if "ener" in tools:

        ener_array = array[:, 0]
        ener = np.mean(ener_array, axis=0)

        for ener_state_val in ener:
            data['ener'].write(f"{ener_state_val}\n")

    #############
    # ener_spac #
    #############

    if "ener_spac" in tools:

        ener_spac_array = array[:, 1, :-1]  # truncate last None value

        for i in range(np.shape(ener_spac_array)[0]):
            for j in range(np.shape(ener_spac_array)[1]):
                data['ener_spac'].write(f"{ener_spac_array[i, j]}\n")

    #######
    # ent #
    #######

    if "ent" in tools:

        ent_array = array[:, 2]
        ent = np.mean(ent_array, axis=0)

        for ent_val in ent:
            data['ent'].write(f"{ent_val}\n")

    ###########
    # ent_mid #
    ###########

    if "ent_mid" in tools:

        ent_mid_array = array[:, 2]
        ent_mid = np.mean(ent_mid_array, axis=0)

        data['ent_mid'].write(f"{_leaf_args['L']}\t{ent_mid}\n")

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('inst_ham')

    my_inst_ham(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
