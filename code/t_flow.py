# --- python imports
import numpy as np
from time import perf_counter
import sys
from joblib import delayed, Parallel
# --- driven_systems imports
import functions.func_ham as fh
import functions.func_args as fa
import functions.func_proc as fp


def bloch_state(cos_theta, phi):

    theta = np.arccos(cos_theta)

    return np.array([np.cos(theta/2), np.exp(1j*phi)*np.sin(theta/2)])


def my_t_flow(path_flag, threads, model, _leaf_args):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = perf_counter()  # start the timer

    leaf = fp.file_name_leaf("t_flow", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("t_flow", path, model, leaf)

    tools = ["ent_t_flow"]
    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    def realization(itr, _t_list, _model, _leaf_args):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        _S_array = np.zeros(_leaf_args['t_samp'], dtype=float)

        H = fh.chosen_hamiltonian(_model, _leaf_args)

        # initial spin product state
        # psi_prod = np.zeros(H.basis.Ns)
        # array_idx_prod = H.basis.index(H.basis.state_to_int('010101'))
        # psi_prod[array_idx_prod] = 1.0
        #
        # print(psi_prod.shape, psi_prod)

        # initial bloch product state
        v = 0
        psi_prod = 1
        for i in range(_leaf_args['L']):
            bloch = bloch_state(np.random.choice([-v, v]), np.random.uniform(0, 2*np.pi))
            psi_prod = np.kron(psi_prod, bloch)

        # print(type(psi_prod), psi_prod.shape, psi_prod)
        psi = H.evolve(psi_prod, 0.0, _t_list)
        # print(psi.shape)
        # print(psi[:, 0].shape)
        for i in range(psi.shape[1]):  # t_samp
            _S_array[i] = float(H.basis.ent_entropy(psi[:, i], sub_sys_A=range(H.basis.L//2))["Sent_A"])

        return _S_array

    ###################################################################################################################

    t_list = np.logspace(_leaf_args['t_min'], _leaf_args['t_max'], _leaf_args['t_samp'])

    array = np.asarray(Parallel(n_jobs=threads)(delayed(realization)(i, t_list, model, leaf_args)
                                                for i in range(leaf_args['dis'])), dtype=object)

    ##############
    # ent_t_flow #
    ##############

    if "ent_t_flow" in tools:

        # print("array.shape = ", array.shape)
        # print(array)
        S = np.mean(array, axis=0)
        # print("S.shape = ", S.shape)

        for i, S_val in enumerate(S):
            data['ent_t_flow'].write(f"{t_list[i]}\t{S[i]}\n")

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('t_flow')

    my_t_flow(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
