# --- python imports
import numpy as np
import time
import sys
from joblib import delayed, Parallel
# --- driven_systems imports
from models.heisenberg import heisenberg
import functions.func_args as fa
import functions.func_proc as fp


def my_inst_ham(path_flag, threads, model, _leaf_args):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = time.time()

    leaf = fp.file_name_leaf("inst_ham", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("inst_ham", path, model, leaf)

    tools = ["ener_spec", "ent_spec"]
    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    def realization(itr, _model, _leaf_args, entropy=False):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")
        if _model == "heisenberg":
            H = heisenberg(_leaf_args['L'], _leaf_args['Nup'], _leaf_args['pauli'],
                           _leaf_args['J'][0], _leaf_args['J'][1], _leaf_args['J'][2], _leaf_args['W'])
        else:
            raise ValueError("model not implemented in inst_ham")
        if entropy:
            _E, psi = H.eigh()
            _S = []
            for i in range(psi.shape[1]):
                _S.append(H.basis.ent_entropy(psi[:, i], sub_sys_A=range(H.basis.L // 2))["Sent_A"])
            return _E, _S
        else:
            _E = H.eigvalsh()
            return _E, None

    ###################################################################################################################

    ent_flag = True if "ent_spec" in tools else False
    array = np.asarray(Parallel(n_jobs=threads)(delayed(realization)(i, model, leaf_args, entropy=ent_flag)
                                                for i in range(leaf_args['dis'])), dtype=object)

    #############
    # ener_spec #
    #############

    if "ener_spec" in tools:

        E_array = array[:, 0]
        E = np.mean(E_array, axis=0)

        for i in E:
            data['ener_spec'].write(f"{i}\n")

    ############
    # ent_spec #
    ############

    if "ent_spec" in tools:

        S_array = array[:, 1]
        S = np.mean(S_array, axis=0)

        for i in S:
            data['ent_spec'].write(f"{i}\n")

    print("Total time taken (seconds) = ", time.time() - t0)


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('inst_ham')

    my_inst_ham(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
