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


def find_eigensystem(_model, _leaf_args):
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
        dt_list = np.array(
            [_leaf_args['T1'] / 2.0, _leaf_args['T0'] / 4.0, _leaf_args['delta'] * _leaf_args['T0'] / 4.0])
        Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, VF=True)
    elif _model == "spin2021_2":
        V, H_1, H_2 = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init, T_init = H_1, 0
        H_list = [V, H_1, H_2]
        dt_list = np.array(
            [_leaf_args['T1'] / 2.0, _leaf_args['T0'] / 4.0, _leaf_args['delta'] * _leaf_args['T0'] / 4.0])
        Floq = Floquet({'H_list': H_list, 'dt_list': dt_list}, VF=True)
    else:
        raise ValueError("model not implemented in eig_U")

    return H_init, T_init, Floq


def my_eig_U(path_flag, threads, model, _leaf_args):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = perf_counter()  # start the timer

    leaf = fp.file_name_leaf("eig_U", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("eig_U", path, model, leaf)

    data = fp.prepare_output_files(["eig_U"], path, model, leaf)

    ###################################################################################################################

    def realization(itr, _model, _leaf_args):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        np.random.seed()

        H_init, T_init, Floq = find_eigensystem(_model, _leaf_args)
        E_init, psi_init = H_init.eigh(time=T_init)
        E, psi = Floq.EF, Floq.VF

        init_array = np.concatenate((E_init[:, None].T, psi_init), axis=0)
        U_array = np.concatenate((E[:, None].T, psi), axis=0)

        return init_array, U_array

    ###################################################################################################################

    array = np.stack(Parallel(n_jobs=threads)(delayed(realization)(i, model, leaf_args)
                                              for i in range(leaf_args['dis'])), axis=0)  # (disorder, tool, state)

    data.create_dataset('eig_init', data=np.around(array[:, 0], 12), compression='gzip', shuffle=True, chunks=True)
    data.create_dataset('eig_U', data=np.around(array[:, 1], 12), compression='gzip', shuffle=True, chunks=True)
    data.close()

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('eig_U')

    my_eig_U(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
