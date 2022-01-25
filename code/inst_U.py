# --- python imports
import numpy as np
from time import perf_counter
import sys
from joblib import delayed, Parallel
from quspin.tools.Floquet import Floquet
import matplotlib.pyplot as plt
# --- driven_systems imports
import functions.func_ham as fh
import functions.func_args as fa
import functions.func_proc as fp


def my_inst_U(path_flag, threads, model, _leaf_args):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = perf_counter()  # start the timer

    leaf = fp.file_name_leaf("inst_U", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("inst_U", path, model, leaf)

    tools = ["q_ener", "q_ener_spac", "floq_struc", "loc_len"]
    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    def realization(itr, _model, _leaf_args, eigenstate=False):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        if model == "ponte2015":
            H = fh.chosen_hamiltonian(_model, _leaf_args)
            H_init, T_init = H, _leaf_args['T1']+_leaf_args['T0']/2
            t_list = np.array([0.0, _leaf_args['T1']]) + np.finfo(float).eps
            dt_list = np.array([_leaf_args['T1'], _leaf_args['T0']])
            Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list},
                           VF=eigenstate, UF=eigenstate, thetaF=eigenstate)
        elif model == "ponte2015_2":
            V, H_0 = fh.chosen_hamiltonian(_model, _leaf_args)
            H_init, T_init = H_0, 0
            H_list = [V, H_0]
            dt_list = np.array([_leaf_args['T1'], _leaf_args['T0']])
            Floq = Floquet({'H_list': H_list, 'dt_list': dt_list}, VF=eigenstate, UF=eigenstate, thetaF=eigenstate)
        elif model == "spin2021":
            H = fh.chosen_hamiltonian(_model, _leaf_args)
            H_init, T_init = H, _leaf_args['T1']/2+_leaf_args['T0']/8
            t_list = np.array([0.0, _leaf_args['T1']/2.0, _leaf_args['T1']/2.0 + _leaf_args['T0']/4.0]) \
                + np.finfo(float).eps
            dt_list = np.array([_leaf_args['T1']/2.0, _leaf_args['T0']/4.0, _leaf_args['delta']*_leaf_args['T0']/4.0])
            Floq = Floquet({'H': H, 't_list': t_list, 'dt_list': dt_list}, VF=eigenstate, UF=eigenstate, thetaF=eigenstate)
        elif model == "spin2021_2":
            V, H_1, H_2 = fh.chosen_hamiltonian(_model, _leaf_args)
            H_init, T_init = H_1, 0
            H_list = [V, H_1, H_2]
            dt_list = np.array([_leaf_args['T1']/2.0, _leaf_args['T0']/4.0, _leaf_args['delta']*_leaf_args['T0']/4.0])
            Floq = Floquet({'H_list': H_list, 'dt_list': dt_list}, VF=eigenstate, UF=eigenstate, thetaF=eigenstate)
        else:
            raise ValueError("model not implemented in inst_U")

        if eigenstate:

            _, alpha = H_init.eigh(time=T_init)

            # # print(Floq.UF)
            # fig = plt.figure()
            # ax = fig.add_subplot(111)
            # mat = ax.matshow(np.abs(Floq.UF), cmap='Greys')
            # ax.set_xlabel("$i$")
            # ax.set_ylabel("$j$")
            # ax.set_title(f"$W={_leaf_args['W']}, T={_leaf_args['T0']}$, $L={_leaf_args['L']}$")
            # cbar = plt.colorbar(mat)
            # cbar.set_label("$|U_{\mathrm{F}, i, j}|$")
            # plt.show()
            #
            # fig = plt.figure()
            # ax = fig.add_subplot(111)
            # ax.plot([np.angle(i)/np.pi for i in Floq.thetaF], '.')
            # ax.set_xlabel("$i$")
            # ax.set_ylabel("$\\theta_{\mathrm{F},i}/\pi$")
            # ax.set_title(f"$W={_leaf_args['W']}, T={_leaf_args['T0']}$, $L={_leaf_args['L']}$")
            # plt.show()

            # quasi-energies
            qE = Floq.EF
            # quasi-energy spacings
            qE_spac = []
            for i in range(len(qE)-1):
                qE_spac.append(qE[i+1]-qE[i])
            qE_spac = np.asarray(qE_spac)

            psi = Floq.VF

            # localization length

            # i_array = unit cell index [0,0,1,1,2,2,...]
            i_array = np.array([val for val in range(H_init.basis.Ns//2) for _ in (0, 1)])

            i_0 = np.zeros(H_init.basis.Ns)
            for j in range(H_init.basis.Ns):
                i_0[j] = i_array.dot(np.abs(psi[:, j])**2)

            loc_len = np.zeros(H_init.basis.Ns)
            for j in range(H_init.basis.Ns):
                loc_len[j] = np.sqrt(np.dot((i_array-i_0[j])**2, np.abs(psi[:, j])**2))

            # A2
            A2 = np.zeros((len(psi), len(alpha)))
            for i_idx in range(len(psi)):
                for alpha_idx in range(len(alpha)):
                    A2[alpha_idx, i_idx] = np.abs(np.dot(psi[:, i_idx], alpha[:, alpha_idx]))**2
            return qE, qE_spac, A2[:, 0], loc_len
        else:
            # quasi-energies
            qE = Floq.EF
            # quasi-energy spacings
            qE_spac = []
            for i in range(len(qE) - 1):
                qE_spac.append(qE[i + 1] - qE[i])
            qE_spac = np.asarray(qE_spac)
            return qE, qE_spac, None, None

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

    ###########
    # loc_len #
    ###########

    if "loc_len" in tools:

        loc_len_array = array[:, 3]
        loc_len = np.mean(loc_len_array, axis=0)

        for i in loc_len:
            data['loc_len'].write(f"{i}\n")

        print(f"average localization length = {np.mean(loc_len)}")

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('inst_U')

    my_inst_U(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
