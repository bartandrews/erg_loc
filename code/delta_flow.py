# --- python imports
import numpy as np
from time import perf_counter
import sys
from joblib import delayed, Parallel
# --- QuSpin imports
from quspin.tools.Floquet import Floquet
# --- driven_systems imports
import functions.func_ham as fh
import functions.func_args as fa
import functions.func_proc as fp


def find_Ns(_model, _leaf_args):

    _leaf_args['delta'] = _leaf_args['delta_min']

    if _model == "ponte2015":
        H = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init = H
    elif _model == "ponte2015_2":
        V, H_0 = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init = H_0
    elif _model == "spin2021":
        H = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init = H
    elif _model == "spin2021_2":
        V, H_1, H_2 = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init = H_1
    else:
        H_init = fh.chosen_hamiltonian(_model, _leaf_args)

    return H_init.Ns


def find_eigensystem(_model, _leaf_args, _delta):
    if _model == "ponte2015":
        H = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init, T_init = H, _leaf_args['T1'] + _leaf_args['T0'] / 2
        t_list = np.array([0.0, _leaf_args['T1']]) + np.finfo(float).eps
        dt_list = np.array([_leaf_args['T1'], _leaf_args['T0']])
        Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, VF=True)
    elif _model == "ponte2015_2":
        V, H_0 = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init, T_init = H_0, 0
        H_list = [V, H_0]
        dt_list = np.array([_leaf_args['T1'], _leaf_args['T0']])
        Floq = Floquet({'H_list': H_list, 'dt_list': dt_list}, VF=True)
    elif _model == "spin2021":
        H = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init, T_init = H, _leaf_args['T1'] / 2 + _leaf_args['T0'] / 8
        t_list = np.array([0.0, _leaf_args['T1'] / 2.0, _leaf_args['T1'] / 2.0 + _leaf_args['T0'] / 4.0]) \
            + np.finfo(float).eps
        dt_list = np.array([_leaf_args['T1'] / 2.0, _leaf_args['T0'] / 4.0, _delta * _leaf_args['T0'] / 4.0])
        Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, VF=True)
    elif _model == "spin2021_2":
        V, H_1, H_2 = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init, T_init = H_1, 0
        H_list = [V, H_1, H_2]
        dt_list = np.array([_leaf_args['T1'] / 2.0, _leaf_args['T0'] / 4.0, _delta * _leaf_args['T0'] / 4.0])
        Floq = Floquet({'H_list': H_list, 'dt_list': dt_list}, VF=True)
    else:
        raise ValueError("model not implemented in delta_flow")

    return H_init, T_init, Floq


def my_delta_flow(path_flag, threads, model, _leaf_args):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = perf_counter()  # start the timer

    leaf = fp.file_name_leaf("delta_flow", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("delta_flow", path, model, leaf)

    # "q_ener_delta_flow", "loc_len_delta_flow", "PR_delta_flow", "ent_delta_flow"
    tools = ["ent_delta_flow"]

    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    def realization(itr, _delta_list, _model, _leaf_args):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        Ns = find_Ns(_model, _leaf_args)

        _q_ener_array = np.zeros((_leaf_args['delta_samp'], Ns))
        _loc_len_array = np.zeros((_leaf_args['delta_samp'], Ns))
        _PR_array = np.zeros((_leaf_args['delta_samp'], Ns))
        _ent_array = np.zeros((_leaf_args['delta_samp'], Ns))

        for i, _delta in enumerate(_delta_list):

            H_init, T_init, Floq = find_eigensystem(_model, _leaf_args, _delta)

            _, alpha = H_init.eigh(time=T_init)
            qE, psi = Floq.EF, Floq.VF

            # --- q_ener_delta_flow
            for k, q_ener_val in enumerate(qE):
                _q_ener_array[i, k] = q_ener_val

            # --- loc_len_delta_flow
            i_array = np.array([k//2 for k in range(H_init.Ns)])
            i_0_array = np.array([i_array.dot(np.abs(psi[:, k])**2) for k in range(H_init.Ns)])
            for k in range(H_init.Ns):
                _loc_len_array[i, k] = np.sqrt(np.dot((i_array-i_0_array[k])**2, np.abs(psi[:, k])**2))

            # --- PR_delta_flow
            A4 = np.zeros((len(alpha), len(psi)))
            for i_idx in range(len(psi)):
                for alpha_idx in range(len(alpha)):
                    A4[alpha_idx, i_idx] = np.abs(np.dot(psi[:, i_idx], alpha[:, alpha_idx]))**4
            for i_idx in range(len(psi)):
                _PR_array[i, i_idx] = (1/H_init.Ns) * (1/np.sum(A4[:, i_idx]))

            # --- ent_delta_flow
            if not any(item in ["loc_len_delta_flow", "PR_delta_flow"] for item in tools):
                for k in range(H_init.Ns):
                    _ent_array[i, k] = H_init.basis.ent_entropy(psi[:, k], sub_sys_A=range(H_init.basis.L//2),
                                                                density=False)["Sent_A"]

        return _q_ener_array, _loc_len_array, _PR_array, _ent_array

    ###################################################################################################################

    delta_list = np.linspace(_leaf_args['delta_min'], _leaf_args['delta_max'], _leaf_args['delta_samp'])
    array = np.stack(Parallel(n_jobs=threads)(delayed(realization)(i, delta_list, model, leaf_args)
                                              for i in range(leaf_args['dis'])), axis=0)
    # (disorder, tool, samp, state)

    #####################
    # q_ener_delta_flow #
    #####################

    if "q_ener_delta_flow" in tools:

        q_ener_array = array[:, 0]
        q_ener = np.mean(q_ener_array, axis=0)

        for i in range(np.shape(q_ener)[0]):
            for j in range(np.shape(q_ener)[1]):
                data['q_ener_delta_flow'].write(f"{delta_list[i]:g}\t{q_ener[i, j]}\n")

    ######################
    # loc_len_delta_flow #
    ######################

    if "loc_len_delta_flow" in tools:

        loc_len_array = array[:, 1]
        loc_len = np.mean(loc_len_array, axis=(0, 2))

        for i, loc_len_samp_val in enumerate(loc_len):
            data['loc_len_delta_flow'].write(f"{delta_list[i]:g}\t{loc_len_samp_val}\n")

    #################
    # PR_delta_flow #
    #################

    if "PR_delta_flow" in tools:

        PR_array = array[:, 2]
        PR = np.mean(PR_array, axis=(0, 2))

        for i, PR_samp_val in enumerate(PR):
            data['PR_delta_flow'].write(f"{delta_list[i]:g}\t{PR_samp_val}\n")

    ##################
    # ent_delta_flow #
    ##################

    if "ent_delta_flow" in tools:

        ent_array = array[:, 3]
        ent = np.mean(ent_array, axis=0)

        for samp in range(np.shape(ent)[0]):
            string = f"{delta_list[samp]:g}"
            for state in range(np.shape(ent)[1]):
                if ent[samp][state] != 0.0:
                    string += f"\t{ent[samp][state]}"
                else:
                    continue
            data['ent_delta_flow'].write(f"{string}\n")

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('delta_flow')

    my_delta_flow(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
