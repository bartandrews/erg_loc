import sys
import os
import ctypes


def file_name_stem(tool, model):
    return f"{tool}_{model}_".replace(" ", "_")


def file_name_leaf(program, model, ham_params):

    L = f"L_{ham_params['L']}_"
    bc = f"{ham_params['bc']}bc_"
    if ham_params['dis'] == 1:
        dis = ""
    else:
        dis = f"dis_{ham_params['dis']}_"
    J = f"J_{ham_params['J'][0]:g}_{ham_params['J'][1]:g}_{ham_params['J'][2]:g}_"
    W = f"W_{ham_params['W']:g}"
    ext = ".dat"

    leaf = f"{L}{bc}{dis}{J}{W}{ext}{ham_params['tag']}"

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
