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

    # "ent_t_flow", "numb_fluc_t_flow", "overlap_t_flow"
    tools = ["overlap_t_flow"]
    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    # noinspection PyUnboundLocalVariable
    def realization(itr, _t_list, _model, _leaf_args):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        _ent_array = np.zeros(_leaf_args['t_samp'])
        _numb_fluc_array = np.zeros(_leaf_args['t_samp'])
        _overlap_array = np.zeros(_leaf_args['t_samp'])

        H = fh.chosen_hamiltonian(_model, _leaf_args)

        initial_state = "custom"

        if initial_state is "bloch":
            v = 1
            psi_prod = 1
            for i in range(_leaf_args['L']):
                bloch = bloch_state(np.random.choice([-v, v]), np.random.uniform(0, 2*np.pi))
                psi_prod = np.kron(psi_prod, bloch)
        elif initial_state is "custom":
            psi_prod = np.zeros(H.Ns)
            idx = H.basis.index(int('10'*12, 2))  # Z2
            # idx = H.basis.index(int('100'*8, 2))  # Z3
            # idx = H.basis.index(int('1000'*6, 2))  # Z4
            psi_prod[idx] = 1
        else:
            raise ValueError("initial state not implemented")

        psi = H.evolve(psi_prod, 0.0, _t_list)

        if "numb_fluc_t_flow" in tools:
            basis = spin_basis_1d(H.basis.L)
            Id_term = [[1, i] for i in range(H.basis.L//2)]
            Sz_term = [[1, i] for i in range(H.basis.L//2)]
            static = [["I", Id_term], ["z", Sz_term]]
            dynamic = []
            total_spin_half_chain = 0.5*hamiltonian(static, dynamic, basis=basis, dtype=np.float64,
                                                    check_symm=False, check_herm=False, check_pcon=False)

        for i in range(_leaf_args['t_samp']):
            if "ent_t_flow" in tools:
                _ent_array[i] = H.basis.ent_entropy(psi[:, i], sub_sys_A=range(H.basis.L//2))["Sent_A"]
            if "numb_fluc_t_flow" in tools:
                _numb_fluc_array[i] = np.real(total_spin_half_chain.quant_fluct(psi[:, i]))
            if "overlap_t_flow" in tools:
                _overlap_array[i] = np.abs(np.dot(psi_prod, psi[:, i]))**2

        return _ent_array, _numb_fluc_array, _overlap_array

    ###################################################################################################################

    t_list = np.linspace(_leaf_args['t_min'], _leaf_args['t_max'], _leaf_args['t_samp'])  # linspace/logspace
    array = np.stack(Parallel(n_jobs=threads)(delayed(realization)(i, t_list, model, leaf_args)
                                              for i in range(leaf_args['dis'])), axis=0)  # (disorder, tool, samp)

    ##############
    # ent_t_flow #
    ##############

    if "ent_t_flow" in tools:

        ent_array = array[:, 0]
        ent = np.mean(ent_array, axis=0)

        for i, ent_val in enumerate(ent):
            data['ent_t_flow'].write(f"{t_list[i]:g}\t{ent_val}\n")

    ####################
    # numb_fluc_t_flow #
    ####################

    if "numb_fluc_t_flow" in tools:

        numb_fluc_array = array[:, 1]
        numb_fluc = np.mean(numb_fluc_array, axis=0)

        for i, numb_fluc_val in enumerate(numb_fluc):
            data['numb_fluc_t_flow'].write(f"{t_list[i]:g}\t{numb_fluc_val}\n")

    ##################
    # overlap_t_flow #
    ##################

    if "overlap_t_flow" in tools:

        overlap_array = array[:, 2]
        overlap = np.mean(overlap_array, axis=0)

        for i, overlap_val in enumerate(overlap):
            data['overlap_t_flow'].write(f"{t_list[i]:g}\t{overlap_val}\n")

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('t_flow')

    my_t_flow(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
