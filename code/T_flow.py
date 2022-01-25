# --- python imports
import numpy as np
from time import perf_counter
import sys
from joblib import delayed, Parallel
from quspin.tools.Floquet import Floquet
# --- driven_systems imports
import functions.func_ham as fh
import functions.func_args as fa
import functions.func_proc as fp


def my_T_flow(path_flag, threads, model, _leaf_args):
    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = perf_counter()  # start the timer

    leaf = fp.file_name_leaf("T_flow", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("T_flow", path, model, leaf)

    tools = ["PR_T_flow"]
    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    def realization(itr, _T_list, _model, _leaf_args):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        _PR = []
        for i, _T in enumerate(_T_list):
            _leaf_args['T1'] = _T
            if model == "ponte2015":
                H = fh.chosen_hamiltonian(_model, _leaf_args)
                H_init, T_init = H, _leaf_args['T1'] + _leaf_args['T0']/2
                t_list = np.array([0.0, _leaf_args['T1']]) + np.finfo(float).eps
                dt_list = np.array([_leaf_args['T1'], _leaf_args['T0']])
                Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, VF=True)
            elif model == "ponte2015_2":
                V, H_0 = fh.chosen_hamiltonian(_model, _leaf_args)
                H_init, T_init = H_0, 0
                H_list = [V, H_0]
                dt_list = np.array([_leaf_args['T1'], _leaf_args['T0']])
                Floq = Floquet({'H_list': H_list, 'dt_list': dt_list}, VF=True)
            elif model == "spin2021":
                H = fh.chosen_hamiltonian(_model, _leaf_args)
                H_init, T_init = H, _leaf_args['T1']/2 + _leaf_args['T0']/8
                t_list = np.array([0.0, _leaf_args['T1']/2.0, _leaf_args['T1']/2.0 + _leaf_args['T0']/4.0]) \
                    + np.finfo(float).eps
                dt_list = np.array([_leaf_args['T1']/2.0, _leaf_args['T0']/4.0,
                                    _leaf_args['delta']*_leaf_args['T0']/4.0])
                Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, VF=True)
            elif model == "spin2021_2":
                V, H_1, H_2 = fh.chosen_hamiltonian(_model, _leaf_args)
                H_init, T_init = H_1, 0
                H_list = [V, H_1, H_2]
                dt_list = np.array([_leaf_args['T1']/2.0, _leaf_args['T0']/4.0,
                                    _leaf_args['delta']*_leaf_args['T0']/4.0])
                Floq = Floquet({'H_list': H_list, 'dt_list': dt_list}, VF=True)
            else:
                raise ValueError("model not implemented in T_flow")

            _, alpha = H_init.eigh(time=T_init)

            # A4
            psi = Floq.VF
            A4 = np.zeros((len(alpha), len(psi)))
            for i_idx in range(len(psi)):
                for alpha_idx in range(len(alpha)):
                    A4[alpha_idx, i_idx] = np.abs(np.dot(psi[:, i_idx], alpha[:, alpha_idx])) ** 4

            _PR_temp = []
            for i_idx in range(len(psi)):
                _PR_temp.append((1/H_init.basis.Ns)*(1/np.sum(A4[:, i_idx])))

            _PR.append(np.mean(_PR_temp))

        return _PR

    ###################################################################################################################

    T_list = np.linspace(_leaf_args['T_min'], _leaf_args['T_max'], _leaf_args['T_samp'])
    array = np.asarray(Parallel(n_jobs=threads)(delayed(realization)(i, T_list, model, leaf_args)
                                                for i in range(leaf_args['dis'])), dtype=object)

    #############
    # PR_T_flow #
    #############

    if "PR_T_flow" in tools:

        PR = np.mean(array, axis=0)

        for i, PR_val in enumerate(PR):
            data['PR_T_flow'].write(f"{T_list[i]:g}\t{PR_val}\n")

    print(f"Total time taken (seconds) = {perf_counter() - t0:.1f}")


if __name__ == "__main__":
    prog_args, stem_args, leaf_args = fa.parse_input_arguments('T_flow')

    my_T_flow(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
