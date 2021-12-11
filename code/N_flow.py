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


def my_N_flow(path_flag, threads, model, _leaf_args):
    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = perf_counter()  # start the timer

    leaf = fp.file_name_leaf("N_flow", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("N_flow", path, model, leaf)

    tools = ["ener_abs_N_flow"]
    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    def realization(itr, _model, _leaf_args):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        if model == "ponte2015":
            H = fh.chosen_hamiltonian(_model, _leaf_args)
            H_init = H
            t_list = np.array([0.0, _leaf_args['T0'] / 2.0, _leaf_args['T0'] / 2.0 + _leaf_args['T1']]) \
                + np.finfo(float).eps
            dt_list = np.array([_leaf_args['T0'] / 2.0, _leaf_args['T1'], _leaf_args['T0'] / 2.0])
            Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, UF=True)
        elif model == "ponte2015_2":
            H_0, V = fh.chosen_hamiltonian(_model, _leaf_args)
            H_init = H_0
            H_list = [H_0, V, H_0]
            dt_list = np.array([_leaf_args['T0'] / 2.0, _leaf_args['T1'], _leaf_args['T0'] / 2.0])
            Floq = Floquet({'H_list': H_list, 'dt_list': dt_list}, UF=True)
        elif model == "spin2021":
            H = fh.chosen_hamiltonian(_model, _leaf_args)
            H_init = H
            t_list = np.array([0.0, _leaf_args['T0'] / 4.0, 3 * _leaf_args['T0'] / 4.0]) + np.finfo(float).eps
            dt_list = np.array([_leaf_args['T0'] / 4.0, _leaf_args['T0'] / 2.0, _leaf_args['T0'] / 4.0])
            Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, UF=True)
        elif model == "spin2021_2":
            V, H_1, H_2 = fh.chosen_hamiltonian(_model, _leaf_args)
            H_init = V
            H_list = [V, H_1, H_2]
            delta = 1
            print("delta = ", delta)
            dt_list = np.array([_leaf_args['T1'] / 2.0, _leaf_args['T0'] / 4.0, delta * _leaf_args['T0'] / 4.0])
            Floq = Floquet({'H_list': H_list, 'dt_list': dt_list}, UF=True)
        else:
            raise ValueError("model not implemented in N_flow")

        # --- energy absorbed under driving
        E_Tinf = H_init.trace(time=0) / H_init.basis.Ns
        E, phi = H_init.eigh(time=0)

        E_0 = np.min(E)
        phi_0 = phi[:, np.argmin(E)]

        UF = Floq.UF

        phi_N = phi_0
        Q_N = np.zeros(_leaf_args['N'])
        for n in range(_leaf_args['N']):

            if n > 0:
                phi_N = UF.dot(phi_N)

            Q_N[n] = (np.real(H_init.matrix_ele(phi_N, phi_N, time=0)) - E_0) / (E_Tinf - E_0)

        return Q_N

    ###################################################################################################################

    array = np.asarray(Parallel(n_jobs=threads)(delayed(realization)(i, model, leaf_args)
                                                for i in range(leaf_args['dis'])), dtype=object)

    ###################
    # ener_abs_N_flow #
    ###################

    if "ener_abs_N_flow" in tools:

        E_abs = np.mean(array, axis=0)

        for i in E_abs:
            data['ener_abs_N_flow'].write(f"{i}\n")

    print(f"Total time taken (seconds) = {perf_counter() - t0:.1f}")


if __name__ == "__main__":
    prog_args, stem_args, leaf_args = fa.parse_input_arguments('N_flow')

    my_N_flow(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
