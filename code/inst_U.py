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


def my_inst_U(path_flag, threads, model, _leaf_args):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = perf_counter()  # start the timer

    leaf = fp.file_name_leaf("inst_U", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("inst_U", path, model, leaf)

    tools = ["q_ener", "q_ener_spac", "floq_struc"]
    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    def realization(itr, _model, _leaf_args, eigenstate=False):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        H = fh.chosen_hamiltonian(_model, _leaf_args)

        if _model == "ponte2015":
            t_list = np.array([0.0, _leaf_args['T0']/2.0, _leaf_args['T0']/2.0 + _leaf_args['T1']]) + np.finfo(float).eps
            dt_list = np.array([_leaf_args['T0']/2.0, _leaf_args['T1'], _leaf_args['T0']/2.0])
        else:
            delta = 0.01  # fraction of T0/4
            t_list = np.array([0.0, _leaf_args['T0']/4.0]) + np.finfo(float).eps
            dt_list = np.array([_leaf_args['T0']/4.0, (_leaf_args['T0']/4.0)*delta])

        if eigenstate:
            _, alpha = H.eigh(time=0)
            Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, VF=True)
            # quasi-energies
            qE = Floq.EF
            # quasi-energy spacings
            qE_spac = []
            for i in range(len(qE)-1):
                qE_spac.append(qE[i+1]-qE[i])
            qE_spac = np.asarray(qE_spac)
            # A2
            psi = Floq.VF
            A2 = np.zeros((len(psi), len(alpha)))
            for i_idx in range(len(psi)):
                for alpha_idx in range(len(alpha)):
                    A2[alpha_idx, i_idx] = np.abs(np.dot(psi[:, i_idx], alpha[:, alpha_idx]))**2
            return qE, qE_spac, A2[:, 0]
        else:
            Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list})
            # quasi-energies
            qE = Floq.EF
            # quasi-energy spacings
            qE_spac = []
            for i in range(len(qE) - 1):
                qE_spac.append(qE[i + 1] - qE[i])
            qE_spac = np.asarray(qE_spac)
            return qE, qE_spac, None

    ###################################################################################################################

    eig_flag = True if "floq_struc" in tools else False
    array = np.asarray(Parallel(n_jobs=threads)(delayed(realization)(i, model, leaf_args, eigenstate=eig_flag)
                                                for i in range(leaf_args['dis'])), dtype=object)

    ##########
    # q_ener #
    ##########

    if "q_ener" in tools:

        E_array = array[:, 0]
        E = np.mean(E_array, axis=0)

        for i in E:
            data['q_ener'].write(f"{i}\n")

    ###############
    # q_ener_spac #
    ###############

    if "q_ener_spac" in tools:

        E_spac_array = array[:, 1]
        E_spac = np.concatenate(E_spac_array).ravel().tolist()

        for i in E_spac:
            data['q_ener_spac'].write(f"{i}\n")

    ##############
    # floq_struc #
    ##############

    if "floq_struc" in tools:

        A2_array = array[:, 2]
        A2 = np.mean(A2_array, axis=0)

        for i in A2:
            data['floq_struc'].write(f"{i}\n")

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('inst_U')

    my_inst_U(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
