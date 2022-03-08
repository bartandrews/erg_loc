# --- python imports
import numpy as np
from time import perf_counter
import sys
from joblib import delayed, Parallel
import matplotlib.pyplot as plt
# --- QuSpin imports
from quspin.tools.Floquet import Floquet
# --- driven_systems imports
import functions.func_ham as fh
import functions.func_args as fa
import functions.func_proc as fp


def find_eigensystem(_model, _leaf_args, _eigenstate):

    if _model == "ponte2015":
        H = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init, T_init = H, _leaf_args['T1'] + _leaf_args['T0'] / 2
        t_list = np.array([0.0, _leaf_args['T1']]) + np.finfo(float).eps
        dt_list = np.array([_leaf_args['T1'], _leaf_args['T0']])
        Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list},
                       VF=_eigenstate, UF=_eigenstate, thetaF=_eigenstate)
    elif _model == "ponte2015_2":
        V, H_0 = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init, T_init = H_0, 0
        H_list = [V, H_0]
        dt_list = np.array([_leaf_args['T1'], _leaf_args['T0']])
        Floq = Floquet({'H_list': H_list, 'dt_list': dt_list},
                       VF=_eigenstate, UF=_eigenstate, thetaF=_eigenstate)
    elif _model == "spin2021":
        H = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init, T_init = H, _leaf_args['T1'] / 2 + _leaf_args['T0'] / 8
        t_list = np.array([0.0, _leaf_args['T1'] / 2.0, _leaf_args['T1'] / 2.0 + _leaf_args['T0'] / 4.0]) \
            + np.finfo(float).eps
        dt_list = np.array(
            [_leaf_args['T1'] / 2.0, _leaf_args['T0'] / 4.0, _leaf_args['delta'] * _leaf_args['T0'] / 4.0])
        Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list},
                       VF=_eigenstate, UF=_eigenstate, thetaF=_eigenstate)
    elif _model == "spin2021_2":
        V, H_1, H_2 = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init, T_init = H_1, 0
        H_list = [V, H_1, H_2]
        dt_list = np.array(
            [_leaf_args['T1'] / 2.0, _leaf_args['T0'] / 4.0, _leaf_args['delta'] * _leaf_args['T0'] / 4.0])
        Floq = Floquet({'H_list': H_list, 'dt_list': dt_list},
                       VF=_eigenstate, UF=_eigenstate, thetaF=_eigenstate)
    else:
        raise ValueError("model not implemented in inst_U")

    return H_init, T_init, Floq


def my_inst_U(path_flag, threads, model, _leaf_args):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = perf_counter()  # start the timer

    leaf = fp.file_name_leaf("inst_U", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("inst_U", path, model, leaf)

    # "q_ener", "q_ener_spac", "floq_struc", "loc_len"
    tools = ["floq_struc"]
    data = fp.prepare_output_files(tools, path, model, leaf)

    plot_unitary = False

    ###################################################################################################################

    def realization(itr, _model, _leaf_args, eigenstate=False):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        H_init, T_init, Floq = find_eigensystem(_model, _leaf_args, eigenstate)

        _q_ener_array = np.zeros(H_init.Ns)
        _q_ener_spac_array = np.zeros(H_init.Ns)
        _floq_struc_array = np.zeros(H_init.Ns)
        _loc_len_array = np.zeros(H_init.Ns)

        if eigenstate:

            qE, psi = Floq.EF, Floq.VF

            if plot_unitary:

                fig = plt.figure()
                ax = fig.add_subplot(111)
                mat = ax.matshow(np.abs(Floq.UF), cmap='Greys')
                ax.set_xlabel("$i$")
                ax.set_ylabel("$j$")
                ax.set_title(f"$W={_leaf_args['W']}, T={_leaf_args['T0']}$, $L={_leaf_args['L']}$")
                cbar = plt.colorbar(mat)
                cbar.set_label("$|U_{\mathrm{F}, i, j}|$")
                plt.show()

                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot([np.angle(i)/np.pi for i in Floq.thetaF], '.')
                ax.set_xlabel("$i$")
                ax.set_ylabel("$\\theta_{\mathrm{F},i}/\pi$")
                ax.set_title(f"$W={_leaf_args['W']}, T={_leaf_args['T0']}$, $L={_leaf_args['L']}$")
                plt.show()

            if "q_ener" in tools:
                for i, q_ener_val in enumerate(qE):
                    _q_ener_array[i] = q_ener_val
            if "q_ener_spac" in tools:
                _q_ener_spac_array[H_init.Ns-1] = None
                for i in range(H_init.Ns-1):
                    _q_ener_spac_array[i] = qE[i+1] - qE[i]
            if "floq_struc" in tools:
                _, alpha = H_init.eigh(time=T_init)
                for i in range(H_init.Ns):
                    _floq_struc_array[i] = np.abs(np.dot(psi[:, 0], alpha[:, i]))**2
            if "loc_len" in tools:
                i_array = np.array([k//2 for k in range(H_init.Ns)])
                i_0_array = np.array([i_array.dot(np.abs(psi[:, k])**2) for k in range(H_init.Ns)])
                for k in range(H_init.Ns):
                    _loc_len_array[k] = np.sqrt(np.dot((i_array-i_0_array[k])**2, np.abs(psi[:, k])**2))

            return _q_ener_array, _q_ener_spac_array, _floq_struc_array, _loc_len_array
        else:
            qE = Floq.EF

            if "q_ener" in tools:
                for i, q_ener_val in enumerate(qE):
                    _q_ener_array[i] = q_ener_val
            if "q_ener_spac" in tools:
                _q_ener_spac_array[H_init.Ns-1] = None
                for i in range(H_init.Ns-1):
                    _q_ener_spac_array[i] = qE[i+1] - qE[i]

            return _q_ener_array, _q_ener_spac_array

    ###################################################################################################################

    if any(item in ["floq_struc", "loc_len"] for item in tools):
        eig_flag = True
    else:
        eig_flag = False

    array = np.stack(Parallel(n_jobs=threads)(delayed(realization)(i, model, leaf_args, eigenstate=eig_flag)
                                              for i in range(leaf_args['dis'])), axis=0)  # (disorder, tool, state)

    ##########
    # q_ener #
    ##########

    if "q_ener" in tools:

        q_ener_array = array[:, 0]
        q_ener = np.mean(q_ener_array, axis=0)

        for q_ener_state_val in q_ener:
            data['q_ener'].write(f"{q_ener_state_val}\n")

    ###############
    # q_ener_spac #
    ###############

    if "q_ener_spac" in tools:

        q_ener_spac_array = array[:, 1, :-1]  # truncate last None value

        for i in range(np.shape(q_ener_spac_array)[0]):
            for j in range(np.shape(q_ener_spac_array)[1]):
                data['q_ener_spac'].write(f"{q_ener_spac_array[i, j]}\n")

    ##############
    # floq_struc #
    ##############

    if "floq_struc" in tools:

        floq_struc_array = array[:, 2]
        floq_struc = np.mean(floq_struc_array, axis=0)

        for floq_struc_state_val in floq_struc:
            data['floq_struc'].write(f"{floq_struc_state_val}\n")

    ###########
    # loc_len #
    ###########

    if "loc_len" in tools:

        loc_len_array = array[:, 3]
        loc_len = np.mean(loc_len_array, axis=0)

        for loc_len_state_val in loc_len:
            data['loc_len'].write(f"{loc_len_state_val}\n")

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('inst_U')

    my_inst_U(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
