import os
import numpy as np

if __name__ == '__main__':

    # loop over configurations
    for L in [10, 12]:
        for delta in np.linspace(0, 1, 11):
            file = open(os.path.join("ent_delta_flow",
                                     f"eig_U_spin2021_L_{L}_obc_dis_100_bat_1_J_1_1_1_T0_1_T1_1_delta_{delta:g}_W_2.sh"),
                        "w+")
            file.write("#!/bin/bash\n"
                       "#$ -l h_rt=3:00:00:00\n"
                       "#$ -l h_data=16G\n"
                       "#$ -l h_vmem=512G\n"
                       "#$ -cwd\n"
                       "#$ -o stdout.$JOB_ID\n"
                       "#$ -j y\n"
                       "#$ -pe shared 32\n"
                       "#$ -l highp\n\n")
            file.write("source /u/home/b/baandr12/.bash_profile\n"
                       "source /u/home/b/baandr12/.bashrc\n\n")
            file.write(f"python code/eig_U.py -mod spin2021 -L {L} -bc o "
                       f"-W 2 -T0 1 -T1 1 -dis 100 -bat 1 -delta {delta}")
            file.close()
