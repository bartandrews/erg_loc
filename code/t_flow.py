# --- python imports
import numpy as np
from time import perf_counter
import sys
from joblib import delayed, Parallel
# --- QuSpin imports
from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian
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

    # "ent_t_flow", "numb_fluc_t_flow"
    tools = ["ent_t_flow"]
    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    def realization(itr, _t_list, _model, _leaf_args):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        _S_array = np.zeros(_leaf_args['t_samp'], dtype=float)
        _numb_fluc_array = np.zeros(_leaf_args['t_samp'], dtype=float)

        H = fh.chosen_hamiltonian(_model, _leaf_args)

        # - initial bloch product state
        v = 1
        psi_prod = 1
        for i in range(_leaf_args['L']):
            bloch = bloch_state(np.random.choice([-v, v]), np.random.uniform(0, 2*np.pi))
            psi_prod = np.kron(psi_prod, bloch)

        psi = H.evolve(psi_prod, 0.0, _t_list)

        basis = spin_basis_1d(H.basis.L)
        Id_term = [[1, i] for i in range(H.basis.L//2)]
        Sz_term = [[1, i] for i in range(H.basis.L//2)]
        static = [["I", Id_term], ["z", Sz_term]]
        dynamic = []
        total_spin_half_chain = 0.5*hamiltonian(static, dynamic, basis=basis, dtype=np.float64,
                                                check_symm=False, check_herm=False, check_pcon=False)

        for i in range(psi.shape[1]):  # t_samp
            _S_array[i] = float(H.basis.ent_entropy(psi[:, i], sub_sys_A=range(H.basis.L//2))["Sent_A"])
            _numb_fluc_array[i] = np.real(total_spin_half_chain.quant_fluct(psi[:, i]))

        return _S_array, _numb_fluc_array

    ###################################################################################################################

    t_list = np.logspace(_leaf_args['t_min'], _leaf_args['t_max'], _leaf_args['t_samp'])

    array = np.asarray(Parallel(n_jobs=threads)(delayed(realization)(i, t_list, model, leaf_args)
                                                for i in range(leaf_args['dis'])), dtype=object)

    ##############
    # ent_t_flow #
    ##############

    if "ent_t_flow" in tools:

        ent_array = array[:, 0]
        S = np.mean(ent_array, axis=0)

        for i, S_val in enumerate(S):
            data['ent_t_flow'].write(f"{t_list[i]:g}\t{S[i]}\n")

    ####################
    # numb_fluc_t_flow #
    ####################

    if "numb_fluc_t_flow" in tools:

        numb_fluc_array = array[:, 1]
        numb_fluc = np.mean(numb_fluc_array, axis=0)

        for i, numb_fluc_val in enumerate(numb_fluc):
            data['numb_fluc_t_flow'].write(f"{t_list[i]:g}\t{numb_fluc[i]}\n")

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('t_flow')

    my_t_flow(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
