import os

if __name__ == '__main__':

    # loop over configurations
    for L in [10, 12, 14, 16]:
        file = open(os.path.join("ent_delta_flow",
                                 f"ent_delta_flow_spin2021_L_{L}_obc_dis_100_J_1_1_1_T0_1_T1_1_delta_0_1_11_W_2.sh"),
                    "w+")
        file.write("#!/bin/bash\n\n")
        file.write("source /u/home/b/baandr12/.bash_profile\n"
                   "source /u/home/b/baandr12/.bashrc\n\n")
        file.write(f"python code/delta_flow.py -mod spin2021 -dis 100 -L {L} -bc o -T0 1 -T1 1 -W 2 "
                   f"-delta_min 0 -delta_max 1 -delta_samp 11")
        file.close()
