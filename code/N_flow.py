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

    if _model == "ising_2":
        H = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init, T_init = H, 0
        Floq = Floquet({'H': H, 'T': _leaf_args['T0']}, UF=True, VF=True)
    elif _model == "ponte2015":
        H = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init, T_init = H, _leaf_args['T1'] + _leaf_args['T0'] / 2
        t_list = np.array([0.0, _leaf_args['T1']]) + np.finfo(float).eps
        dt_list = np.array([_leaf_args['T1'], _leaf_args['T0']])
        Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, UF=True, VF=True)
    elif _model == "ponte2015_2":
        V, H_0 = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init, T_init = H_0, 0
        H_list = [V, H_0]
        dt_list = np.array([_leaf_args['T1'], _leaf_args['T0']])
        Floq = Floquet({'H_list': H_list, 'dt_list': dt_list}, UF=True, VF=True)
    elif _model == "spin2021":
        H = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init, T_init = H, _leaf_args['T1'] / 2 + _leaf_args['T0'] / 8
        t_list = np.array([0.0, _leaf_args['T1'] / 2.0, _leaf_args['T1'] / 2.0 + _leaf_args['T0'] / 4.0]) \
            + np.finfo(float).eps
        dt_list = np.array(
            [_leaf_args['T1'] / 2.0, _leaf_args['T0'] / 4.0, _leaf_args['delta'] * _leaf_args['T0'] / 4.0])
        Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, UF=True, VF=True)
    elif _model == "spin2021_2":
        V, H_1, H_2 = fh.chosen_hamiltonian(_model, _leaf_args)
        H_init, T_init = H_1, 0
        H_list = [V, H_1, H_2]
        dt_list = np.array(
            [_leaf_args['T1'] / 2.0, _leaf_args['T0'] / 4.0, _leaf_args['delta'] * _leaf_args['T0'] / 4.0])
        Floq = Floquet({'H_list': H_list, 'dt_list': dt_list}, UF=True, VF=True)
    else:
        raise ValueError("model not implemented in N_flow")

    return H_init, T_init, Floq


def bloch_state(cos_theta, phi):

    theta = np.arccos(cos_theta)

    return np.array([np.cos(theta/2), np.exp(1j*phi)*np.sin(theta/2)])


def my_N_flow(path_flag, threads, model, _leaf_args):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = perf_counter()  # start the timer

    leaf = fp.file_name_leaf("N_flow", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("N_flow", path, model, leaf)

    # "ener_abs_N_flow", "ent_N_flow", "ent_info_N_flow"
    tools = ["ent_info_N_flow"]
    if "ent_info_N_flow" in tools and len(tools) != 1:
        raise ValueError("The tool ent_info_N_flow can only be used in isolation.")

    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    def realization(itr, _model, _leaf_args):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        H_init, T_init, Floq = find_eigensystem(_model, _leaf_args)

        _ener_abs_array = np.zeros(_leaf_args['N'])
        _ent_array = np.zeros(_leaf_args['N'])
        _ent_info_array = np.zeros(_leaf_args['N'])

        E, phi = H_init.eigh(time=T_init)
        UF = Floq.UF

        if any(item in ["ener_abs_N_flow", "ent_N_flow"] for item in tools):
            E_Tinf = H_init.trace(time=T_init) / H_init.Ns
            E_0 = np.min(E)
            phi_0 = phi[:, np.argmin(E)]

            if "ent_N_flow" in tools:
                v = 0
                psi_prod = 1
                for i in range(_leaf_args['L']):
                    bloch = bloch_state(np.random.choice([-v, v]), np.random.uniform(0, 2*np.pi))
                    psi_prod = np.kron(psi_prod, bloch)
                psi_prod_N = psi_prod
            else:  # ener_abs_N_flow
                psi_prod_N = phi_0

            phi_N = phi_0
            for n in range(_leaf_args['N']):
                if n > 0:
                    phi_N = UF.dot(phi_N)
                    psi_prod_N = UF.dot(psi_prod_N)
                _ener_abs_array[n] = (np.real(H_init.matrix_ele(phi_N, phi_N, time=T_init))-E_0)/(E_Tinf-E_0)
                _ent_array[n] = float(H_init.basis.ent_entropy(psi_prod_N,
                                                               sub_sys_A=range(H_init.basis.L//2))["Sent_A"])

        if "ent_info_N_flow" in tools:
            UF_new = UF
            VF = Floq.VF
            for N in range(1, _leaf_args['N']+1):

                if N > 1:
                    UF_new = UF @ UF_new
                    _, VF = np.linalg.eigh(UF_new)

                c = np.zeros((H_init.Ns, H_init.Ns), dtype=complex)
                for n in range(H_init.Ns):
                    for m in range(H_init.Ns):
                        c[n, m] = np.dot(VF[:, n], phi[:, m])

                S = np.zeros(H_init.Ns)
                for n in range(H_init.Ns):
                    c2 = np.square(np.abs(c[n, :]))
                    S[n] = -np.sum(c2*np.log(c2))
                _ent_info_array[N-1] = np.mean(S)/(np.log(0.48*H_init.Ns))

        return _ener_abs_array, _ent_array, _ent_info_array

    ###################################################################################################################

    array = np.stack(Parallel(n_jobs=threads)(delayed(realization)(i, model, leaf_args)
                                              for i in range(leaf_args['dis'])), axis=0)  # (disorder, tool, samp)

    ###################
    # ener_abs_N_flow #
    ###################

    if "ener_abs_N_flow" in tools:

        ener_abs_array = array[:, 0]
        ener_abs = np.mean(ener_abs_array, axis=0)

        for ener_abs_samp_val in ener_abs:
            data['ener_abs_N_flow'].write(f"{ener_abs_samp_val}\n")

    ##############
    # ent_N_flow #
    ##############

    if "ent_N_flow" in tools:

        ent_array = array[:, 1]
        ent = np.mean(ent_array, axis=0)

        for ent_samp_val in ent:
            data['ent_N_flow'].write(f"{ent_samp_val}\n")

    ###################
    # ent_info_N_flow #
    ###################

    if "ent_info_N_flow" in tools:

        ent_info_array = array[:, 2]
        ent_info = np.mean(ent_info_array, axis=0)

        for ent_info_samp_val in ent_info:
            data['ent_info_N_flow'].write(f"{ent_info_samp_val}\n")

    print(f"Total time taken (seconds) = {perf_counter() - t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('N_flow')

    my_N_flow(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
