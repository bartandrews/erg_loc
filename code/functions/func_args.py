# --- python imports
import argparse


def parse_input_arguments(program):

    parser = argparse.ArgumentParser(prog=program)
    prog = parser.add_argument_group("program sub-arguments")
    stem = parser.add_argument_group("stem sub-arguments")
    leaf = parser.add_argument_group("leaf sub-arguments")

    models = ["heisenberg"]

    # program sub-arguments
    prog.add_argument("-path", default=False, action='store_true', help="use a custom path")
    prog.add_argument("-thr", "--threads", type=int, default=-1, help="number of threads")

    # stem sub-arguments (model specific)
    stem.add_argument("-mod", "--model", type=str, default="heisenberg", choices=models, required=True,
                      help="name of model")

    # leaf sub-arguments (run specific)
    leaf.add_argument("-L", type=int, default=8, required=True, help="length of chain")
    leaf.add_argument("-Nup", type=int, default=None, help="number of up spins")
    leaf.add_argument("-pauli", type=int, default=1, choices=[0, 1, -1], help="type of spin operator")
    boundary_conditions = ["o", "p"]
    leaf.add_argument("-bc", type=str, default="o", choices=boundary_conditions, required=True,
                      help="boundary conditions")
    leaf.add_argument("-dis", type=int, default=1, help="number of disorders")
    leaf.add_argument("-J", nargs=3, type=float, default=[1, 1, 1], help="coupling strength (heisenberg model)")
    leaf.add_argument("-W", type=float, default=1, help="disorder strength")
    leaf.add_argument("-tag", type=str, default="", help="tag to append to filename")

    args = vars(parser.parse_args())

    prog_args, stem_args, leaf_args = dict(), dict(), dict()
    for prog_key in ['path', 'threads']:
        prog_args.update({prog_key: args.pop(prog_key, None)})
    for stem_key in ['model']:
        stem_args.update({stem_key: args.pop(stem_key, None)})
    leaf_args = args

    return prog_args, stem_args, leaf_args
