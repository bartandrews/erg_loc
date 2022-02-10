# --- python imports
import numpy as np
from time import perf_counter
import sys
from joblib import delayed, Parallel
# --- driven_systems imports
import functions.func_ham as fh
import functions.func_args as fa
import functions.func_proc as fp
from models.heisenberg import heisenberg


def my_L_flow(path_flag, threads, model, _leaf_args):

    path = "/data/baandr" if path_flag else ""  # specify the custom path
    t0 = perf_counter()  # start the timer

    leaf = fp.file_name_leaf("L_flow", model, _leaf_args)
    sys.stdout = sys.stderr = fp.Logger("L_flow", path, model, leaf)

    # "ent_mid_L_flow"
    tools = ["ent_mid_L_flow"]
    data = fp.prepare_output_files(tools, path, model, leaf)

    ###################################################################################################################

    def realization(itr, _L_list, _model, _leaf_args):
        print(f"Iteration {itr + 1} of {_leaf_args['dis']}")

        _S_array = np.zeros(_leaf_args['L_samp'], dtype=object)
        _S_page_array = np.zeros(_leaf_args['L_samp'], dtype=object)

        _Nup_list = [None]*len(_L_list)
        if _leaf_args['Nup_min'] is not None:
            _Nup_list = list(map(int, np.linspace(_leaf_args['Nup_min'], _leaf_args['Nup_max'], _leaf_args['L_samp'])))

        for i, (_L, _Nup) in enumerate(zip(_L_list, _Nup_list)):
            _leaf_args['L'] = _L
            _leaf_args['Nup'] = _Nup
            H = fh.chosen_hamiltonian(_model, _leaf_args)

            d = heisenberg(_L//2, None, _leaf_args['pauli'],
                           _leaf_args['J'][0], _leaf_args['J'][1], _leaf_args['J'][2],
                           _leaf_args['W']).Ns
            _S_page_array[i] = 0
            for k in range(d+1, d*d+1):
                _S_page_array[i] += 1/k
            _S_page_array[i] -= (d-1)/(2*d)

            Emin, Emax = H.eigsh(k=2, which="BE", maxiter=1E4, return_eigenvectors=False)
            Etarget = Emin + 0.5*(Emax - Emin)
            _, psi = H.eigsh(k=1, sigma=Etarget, maxiter=1E4)
            _S_array[i] = float(H.basis.ent_entropy(psi, sub_sys_A=range(H.basis.L//2), density=False)["Sent_A"])

        return _S_array, _S_page_array

    ###################################################################################################################

    L_list = list(map(int, np.linspace(_leaf_args['L_min'], _leaf_args['L_max'], _leaf_args['L_samp'])))

    array = np.asarray(Parallel(n_jobs=threads)(delayed(realization)(i, L_list, model, leaf_args)
                                                for i in range(leaf_args['dis'])), dtype=object)

    ##################
    # ent_mid_L_flow #
    ##################

    if "ent_mid_L_flow" in tools:

        ent_array = array[:, 0]
        S = np.mean(ent_array, axis=0)

        ent_page_array = array[:, 1]
        S_page = np.mean(ent_page_array, axis=0)

        for i, S_val in enumerate(S):
            data['ent_mid_L_flow'].write(f"{L_list[i]}\t{S_val}\t{S_page[i]}\n")

    print(f"Total time taken (seconds) = {perf_counter()-t0:.1f}")


if __name__ == "__main__":

    prog_args, stem_args, leaf_args = fa.parse_input_arguments('L_flow')

    my_L_flow(prog_args['path'], prog_args['threads'], stem_args['model'], leaf_args)
