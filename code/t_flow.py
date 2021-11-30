# --- python imports
import numpy as np
from time import perf_counter
import sys
from joblib import delayed, Parallel
# --- driven_systems imports
import functions.func_ham as fh
import functions.func_args as fa
import functions.func_proc as fp


def my_t_flow(path_flag, threads, model, _leaf_args):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = perf_counter()  # start the timer

    leaf = fp.file_name_leaf("t_flow", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("t_flow", path, model, leaf)

    tools = ["ent_t_flow"]
    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    def realization(itr, _model, _leaf_args):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        _S_array = np.zeros(_leaf_args['t_samp'], dtype=float)

        t_list = np.linspace(_leaf_args['t_min'], _leaf_args['t_max'], _leaf_args['t_samp'])
        H = fh.chosen_hamiltonian(_model, _leaf_args)
        _E, psi = H.eigh()
        # initial state
        psi_mid = psi[:, np.argsort(_E)[len(_E)//2]]
        # time evolution
        psi = H.evolve(psi_mid, 0.0, t_list)
        for i in range(psi.shape[1]):
            _S_array[i] = float(H.basis.ent_entropy(psi[:, i], sub_sys_A=range(H.basis.L//2))["Sent_A"])

        return _S_array

    ###################################################################################################################

    array = np.asarray(Parallel(n_jobs=threads)(delayed(realization)(i, model, leaf_args)
                                                for i in range(leaf_args['dis'])), dtype=object)

    t_list = np.linspace(_leaf_args['t_min'], _leaf_args['t_max'], _leaf_args['t_samp'])

    ##############
    # ent_t_flow #
    ##############

    if "ent_t_flow" in tools:

        S = np.mean(array, axis=0)

        for i, S_val in enumerate(S):
            data['ent_t_flow'].write(f"{t_list[i]}\t{S[i]}\n")

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('t_flow')

    my_t_flow(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
