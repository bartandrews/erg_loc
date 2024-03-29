# --- python imports
import argparse


def parse_input_arguments(program):

    parser = argparse.ArgumentParser(prog=program)
    prog = parser.add_argument_group("program sub-arguments")
    stem = parser.add_argument_group("stem sub-arguments")
    leaf = parser.add_argument_group("leaf sub-arguments")

    models = ["heisenberg", "ising_2", "ponte2015", "ponte2015_2", "spin2021", "spin2021_2", "pxp"]

    # program sub-arguments

    prog.add_argument("-path", default=False, action='store_true', help="use a custom path")
    prog.add_argument("-thr", "--threads", type=int, default=-1, help="number of threads")

    # stem sub-arguments (model specific)

    stem.add_argument("-mod", "--model", type=str, default="heisenberg", choices=models, required=True,
                      help="name of model")

    # leaf sub-arguments (run specific)

    if program == "L_flow":
        leaf.add_argument("-L_min", type=int, default=8, required=True, help="minimum chain length")
        leaf.add_argument("-L_max", type=int, default=12, required=True, help="maximum chain length")
        leaf.add_argument("-L_samp", type=int, default=3, required=True, help="number of chain length samples")
    else:
        leaf.add_argument("-L", type=int, default=8, required=True, help="length of chain")

    if program == "L_flow":
        leaf.add_argument("-Nup_min", type=int, default=None, help="minimum number of up spins")
        leaf.add_argument("-Nup_max", type=int, default=None, help="maximum number of up spins")
    else:
        leaf.add_argument("-Nup", type=int, default=None, help="number of up spins")

    leaf.add_argument("-pauli", type=int, default=1, choices=[0, 1, -1], help="type of spin operator")
    boundary_conditions = ["o", "p"]
    leaf.add_argument("-bc", type=str, default="o", choices=boundary_conditions, required=True,
                      help="boundary conditions")
    leaf.add_argument("-dis", type=int, default=1, help="number of disorders")
    if program == "eig_U":
        leaf.add_argument("-bat", "--batch", type=int, default=1, help="batch number")

    if program == "t_flow":
        leaf.add_argument("-t_min", type=float, default=0, required=True, help="start time")
        leaf.add_argument("-t_max", type=float, default=10, required=True, help="end time")
        leaf.add_argument("-t_samp", type=int, default=11, required=True, help="number of time samples")

    leaf.add_argument("-J", nargs=3, type=float, default=[1, 1, 1], help="coupling strength")

    leaf.add_argument("-h0", type=float, default=None, help="middle kick strength")
    leaf.add_argument("-T0", type=float, default=None, help="drive interval 0")
    leaf.add_argument("-T1", type=float, default=None, help="drive interval 1")

    if program == "N_flow":
        leaf.add_argument("-N", type=int, default=31, required=True, help="number of Floquet cycles")

    if program == "T_flow":
        leaf.add_argument("-T_min", type=float, default=0, required=True, help="minimum T")
        leaf.add_argument("-T_max", type=float, default=3, required=True, help="maximum T")
        leaf.add_argument("-T_samp", type=int, default=11, required=True, help="number of T samples")

    if program == "W_flow":
        leaf.add_argument("-W_min", type=float, default=0, required=True, help="minimum disorder strength")
        leaf.add_argument("-W_max", type=float, default=10, required=True, help="maximum disorder strength")
        leaf.add_argument("-W_samp", type=int, default=11, required=True, help="number of disorder strength samples")
    else:
        leaf.add_argument("-W", type=float, default=0, help="disorder strength")

    if program == "delta_flow":
        leaf.add_argument("-delta_min", type=float, default=0, required=True, help="minimum delta")
        leaf.add_argument("-delta_max", type=float, default=1, required=True, help="maximum delta")
        leaf.add_argument("-delta_samp", type=int, default=11, required=True, help="number of delta samples")
    else:
        leaf.add_argument("-delta", type=float, default=1, help="fraction of drive (see model for details)")

    leaf.add_argument("-tag", type=str, default="", help="tag to append to filename")

    args = vars(parser.parse_args())

    prog_args, stem_args, leaf_args = dict(), dict(), dict()
    for prog_key in ['path', 'threads']:
        prog_args.update({prog_key: args.pop(prog_key, None)})
    for stem_key in ['model']:
        stem_args.update({stem_key: args.pop(stem_key, None)})
    leaf_args = args

    return prog_args, stem_args, leaf_args
