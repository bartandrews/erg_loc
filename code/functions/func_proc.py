# --- python imports
import sys
import os


def file_name_stem(tool, model):
    return f"{tool}_{model}_".replace(" ", "_")


def file_name_leaf(program, model, ham_params):

    if program == "L_flow":
        L = f"L_{ham_params['L_min']:g}_{ham_params['L_max']:g}_{ham_params['L_samp']}_"
        Nup = f"Nup_{ham_params['Nup_min']}_{ham_params['Nup_max']}_{ham_params['L_samp']}_" if ham_params['Nup_min'] is not None else ""
    else:
        L = f"L_{ham_params['L']:g}_"
        Nup = f"Nup_{ham_params['Nup']}_" if ham_params['Nup'] is not None else ""

    if ham_params['pauli'] != 1:
        pauli = f"pauli_{ham_params['pauli']}_"
    else:
        pauli = ""
    bc = f"{ham_params['bc']}bc_"
    if ham_params['dis'] == 1:
        dis = ""
    else:
        dis = f"dis_{ham_params['dis']:g}_"

    if program == "t_flow":
        t = f"t_{ham_params['t_min']:g}_{ham_params['t_max']:g}_{ham_params['t_samp']}_"
    else:
        t = ""

    J = f"J_{ham_params['J'][0]:g}_{ham_params['J'][1]:g}_{ham_params['J'][2]:g}_"

    if "h0" in ham_params and ham_params['h0'] is not None:
        h0 = f"h0_{ham_params['h0']:g}_"
    else:
        h0 = ""

    if "T0" in ham_params and ham_params['T0'] is not None:
        T0 = f"T0_{ham_params['T0']:g}_"
    else:
        T0 = ""

    if "T1" in ham_params and ham_params['T1'] is not None:
        T1 = f"T1_{ham_params['T1']:g}_"
    else:
        T1 = ""

    if "N" in ham_params:
        N = f"N_{ham_params['N']:g}_"
    else:
        N = ""

    if program == "T_flow":
        T = f"T_{ham_params['T_min']:g}_{ham_params['T_max']:g}_{ham_params['T_samp']}_"
    else:
        T = ""

    if program == "W_flow":
        W = f"W_{ham_params['W_min']:g}_{ham_params['W_max']:g}_{ham_params['W_samp']}"
    else:
        W = f"W_{ham_params['W']:g}"

    ext = ".dat"

    leaf = f"{L}{Nup}{pauli}{bc}{dis}{t}{J}{h0}{T0}{T1}{N}{T}{W}{ext}{ham_params['tag']}"

    return leaf


class Logger(object):
    def __init__(self, program, path, model, leaf):
        self.terminal = sys.stdout or sys.stderr
        os.makedirs(os.path.join(path, "logs", f"{program}", f"{model}", ""), exist_ok=True)
        self.log = open(os.path.join(path, "logs", f"{program}", f"{model}",
                                     f"log_{program}_{model}_{leaf}"), 'w', buffering=1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass


def prepare_output_files(tools, path, model, leaf):

    stem, file, data = [dict()]*3
    for tool in tools:
        stem.update({tool: file_name_stem(tool, model)})
        os.makedirs(os.path.join(path, "data", f"{tool}", f"{model}", ""), exist_ok=True)
        file.update({tool: os.path.join(path, "data", f"{tool}", f"{model}", stem[tool] + leaf)})
        open(file[tool], "w")
        data[tool] = open(file[tool], "a", buffering=1)

    return data
