import numpy as np
import csv
import os

proj_root = '/home/bart/PycharmProjects/erg_loc'

input_file_1 = "data/q_ener/spin2021/q_ener_spin2021_L_20_Nup_1_obc_J_1_1_1_T0_1_T1_1_delta_0.1_W_0.1.dat"
input_file_2 = "data/q_ener/spin2021/q_ener_spin2021_L_20_Nup_2_obc_J_1_1_1_T0_1_T1_1_delta_0.1_W_0.1.dat"

with open(os.path.join(proj_root, input_file_1), 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    E1 = []
    for i, row in enumerate(plots):
        E1.append(float(row[0]))

print("len(E1) = ", len(E1))

# delta = 0
# Omega = 2*np.pi / (0.5+0.25+delta*0.25)

with open(os.path.join(proj_root, input_file_2), 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter='\t')
    E2 = []
    for i, row in enumerate(plots):
        E2.append(float(row[0]))

print("len(E2) = ", len(E2))

indices_to_delete = []
for i, E2_val in enumerate(E2):
    for j1, E1_val1 in enumerate(E1):
        if np.abs(np.abs(E2_val) - np.abs(E1_val1)) < 1e-10:
            if i not in indices_to_delete:
                indices_to_delete.append(i)
        for j2, E1_val2 in enumerate(E1):
            # pos_sup = E1_val1 + E1_val2
            # if pos_sup > Omega:
            #     print("positive exceeded")
            #     pos_sup = - Omega + (pos_sup - Omega)
            # elif pos_sup < -Omega:
            #     pos_sup = Omega - (Omega - pos_sup)
            # neg_sup = E1_val1 - E1_val2
            # if neg_sup < -Omega:
            #     neg_sup = Omega - (Omega - neg_sup)
            # elif neg_sup > Omega:
            #     neg_sup = -Omega + (neg_sup - Omega)
            if np.abs(np.abs(E2_val) - np.abs(E1_val1 + E1_val2)) < 1e-10 \
                    or np.abs(np.abs(E2_val) - np.abs(E1_val1 - E1_val2)) < 1e-10:
                if i not in indices_to_delete:
                    indices_to_delete.append(i)

print(len(indices_to_delete), indices_to_delete)
