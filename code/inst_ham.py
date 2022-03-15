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


def my_inst_ham(path_flag, threads, model, _leaf_args):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = perf_counter()  # start the timer

    leaf = fp.file_name_leaf("inst_ham", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("inst_ham", path, model, leaf)

    # "ener", "ener_spac", "ent", "ent_mid", "overlap", "exp_val", "spec_func"
    tools = ["spec_func"]
    if "ent_mid" in tools and len(tools) != 1:
        raise ValueError("The tool ent_mid can only be used in isolation.")
    if "spec_func" in tools and len(tools) != 1:
        raise ValueError("The tool spec_func can only be used in isolation.")

    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    def realization(itr, _model, _leaf_args, entropy=0):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        H = fh.chosen_hamiltonian(_model, _leaf_args)

        if _model == "pxp":  # user_basis
            L = H.basis.N
        else:
            L = H.basis.L

        _ener_array = np.zeros(H.Ns)
        _ener_spac_array = np.zeros(H.Ns)
        _ent_array = np.zeros(H.Ns)
        _ent_mid = 0
        _overlap_array = np.zeros(H.Ns)
        _exp_val_array = np.zeros(H.Ns)

        if entropy == 2:  # ent
            E, psi = H.eigh()

            if "ener" in tools:
                for i, ener_val in enumerate(E):
                    _ener_array[i] = ener_val
            if "ener_spac" in tools:
                _ener_spac_array[H.Ns-1] = None
                for i in range(H.Ns-1):
                    _ener_spac_array[i] = E[i+1]-E[i]
            if "ent" in tools:
                for i in range(H.Ns):
                    _ent_array[i] = H.basis.ent_entropy(psi[:, i],
                                                        sub_sys_A=range(L//2), density=False)["Sent_A"]
            if "overlap" in tools:
                Z2_state = np.zeros(H.Ns)
                Z2_state[0] = 1
                for i in range(H.Ns):
                    _overlap_array[i] = np.abs(np.dot(Z2_state, psi[:, i]))**2
            if "exp_val" in tools:
                z_term = [[1, i] for i in range(L)]
                static = [["z", z_term]]
                dynamic = []
                Z1 = (1/L)*hamiltonian(static, dynamic, basis=H.basis, dtype=np.float64,
                                       check_symm=False, check_herm=False, check_pcon=False)
                for i in range(H.Ns):
                    _exp_val_array[i] = Z1.expt_value(psi[:, i])
            if "spec_func" in tools:
                _omega_array = np.zeros((2*H.Ns//3, 2*H.Ns//3))
                _spec_func_array = np.zeros((2*H.Ns//3, 2*H.Ns//3))
                z_term = [[1, i] for i in range(L)]
                static = [["z", z_term]]
                dynamic = []
                Z1 = (1/L) * hamiltonian(static, dynamic, basis=H.basis, dtype=np.float64,
                                         check_symm=False, check_herm=False, check_pcon=False)
                if 5*len(E)//6 - len(E)//6 > 2*H.Ns//3:
                    upper_limit = 5*len(E)//6 - 1
                else:
                    upper_limit = 5*len(E)//6
                for a_idx, alpha in enumerate(range(len(E)//6, upper_limit)):
                    for b_idx, beta in enumerate(range(len(E)//6, upper_limit)):
                        _omega_array[a_idx, b_idx] = E[alpha] - E[beta]
                        # e^S(E) prefactor is not yet implemented
                        _spec_func_array[a_idx, b_idx] = np.abs(np.dot(psi[:, beta], Z1.dot(psi[:, alpha])))**2

                _omega_array = _omega_array.ravel()
                _spec_func_array = _spec_func_array.ravel()

                return _omega_array, _spec_func_array

            return _ener_array, _ener_spac_array, _ent_array, _overlap_array, _exp_val_array
        elif entropy == 1:  # ent_mid
            _, psi = H.eigsh(k=1, sigma=0.5, maxiter=1E4)

            if "ent_mid" in tools:
                _ent_mid = H.basis.ent_entropy(psi, sub_sys_A=range(L//2), density=False)["Sent_A"]

            return None, None, _ent_mid
        elif entropy == 0:  # no ent
            E = H.eigvalsh()

            if "ener" in tools:
                for i, ener_val in enumerate(E):
                    _ener_array[i] = ener_val
            if "ener_spac" in tools:
                _ener_spac_array[H.Ns - 1] = None
                for i in range(H.Ns - 1):
                    _ener_spac_array = E[i + 1] - E[i]

            return _ener_array, _ener_spac_array, None
        else:
            raise ValueError("ent_flag not defined in realization")

    ###################################################################################################################

    ent_flag = 0
    if any(item in ["ent", "ent_mid", "overlap", "exp_val", "spec_func"] for item in tools):
        if "ent_mid" in tools:
            ent_flag = 1
        else:
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

        for ent_state_val in ent:
            data['ent'].write(f"{ent_state_val}\n")

    ###########
    # ent_mid #
    ###########

    if "ent_mid" in tools:

        ent_mid_array = array[:, 2]
        ent_mid = np.mean(ent_mid_array, axis=0)

        data['ent_mid'].write(f"{_leaf_args['L']}\t{ent_mid}\n")

    ###########
    # overlap #
    ###########

    if "overlap" in tools:

        overlap_array = array[:, 3]
        overlap = np.mean(overlap_array, axis=0)

        for overlap_val in overlap:
            data['overlap'].write(f"{overlap_val}\n")

    ###########
    # exp_val #
    ###########

    if "exp_val" in tools:

        exp_val_array = array[:, 4]
        exp_val = np.mean(exp_val_array, axis=0)

        for exp_val_val in exp_val:
            data['exp_val'].write(f"{exp_val_val}\n")

    #############
    # spec_func #
    #############

    if "spec_func" in tools:

        omega_array = array[:, 0]
        omega = np.mean(omega_array, axis=0)
        spec_func_array = array[:, 1]
        spec_func = np.mean(spec_func_array, axis=0)

        for i, omega_val in enumerate(omega):
            data['spec_func'].write(f"{omega_val}\t{spec_func[i]}\n")

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('inst_ham')

    my_inst_ham(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
