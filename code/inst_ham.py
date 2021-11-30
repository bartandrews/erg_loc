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

    tools = ["ener", "ener_spac", "ent"]
    if "ent_mid" not in tools:
        data = fp.prepare_output_files(tools, path, model, leaf)
    else:
        data = fp.prepare_output_files(["ent_mid"], path, model, leaf)

    ###################################################################################################################

    def realization(itr, _model, _leaf_args, entropy=0):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        H = fh.chosen_hamiltonian(_model, _leaf_args)

        if entropy == 2:  # ent
            _E, psi = H.eigh()
            _E_spac = []
            for i in range(len(_E)-1):
                _E_spac.append(_E[i+1]-_E[i])
            _E_spac = np.asarray(_E_spac)
            _S = []
            for i in range(psi.shape[1]):
                _S.append(H.basis.ent_entropy(psi[:, i], sub_sys_A=range(H.basis.L // 2))["Sent_A"])
            return _E, _E_spac, _S
        elif entropy == 1:  # ent_mid
            _, psi = H.eigsh(k=1, sigma=0.5, maxiter=1E4)
            _S = H.basis.ent_entropy(psi, sub_sys_A=range(H.basis.L//2))["Sent_A"]
            return None, None, _S
        elif entropy == 0:  # no ent
            _E = H.eigvalsh()
            _E_spac = []
            for i in range(len(_E)-1):
                _E_spac.append(_E[i+1]-_E[i])
            _E_spac = np.asarray(_E_spac)
            return _E, _E_spac, None
        else:
            raise ValueError("ent_flag not defined in realization")

    ###################################################################################################################

    ent_flag = 0
    if any("ent" in i for i in tools):
        if "ent_mid" in tools:
            ent_flag = 1
        elif "ent" in tools:
            ent_flag = 2

    array = np.asarray(Parallel(n_jobs=threads)(delayed(realization)(i, model, leaf_args, entropy=ent_flag)
                                                for i in range(leaf_args['dis'])), dtype=object)

    ########
    # ener #
    ########

    if "ener" in tools and "ent_mid" not in tools:

        E_array = array[:, 0]
        E = np.mean(E_array, axis=0)

        for i in E:
            data['ener'].write(f"{i}\n")

    #############
    # ener_spac #
    #############

    if "ener_spac" in tools and "ent_mid" not in tools:

        E_spac_array = array[:, 1]
        E_spac = np.concatenate(E_spac_array).ravel().tolist()

        for i in E_spac:
            data['ener_spac'].write(f"{i}\n")

    #######
    # ent #
    #######

    if "ent" in tools and "ent_mid" not in tools:

        S_array = array[:, 2]
        S = np.mean(S_array, axis=0)

        for i in S:
            data['ent'].write(f"{i}\n")

    ###########
    # ent_mid #
    ###########

    if "ent_mid" in tools:

        S_array = array[:, 2]
        S = np.mean(S_array, axis=0)

        data['ent_mid'].write(f"{_leaf_args['L']}\t{S}\n")

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('inst_ham')

    my_inst_ham(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
