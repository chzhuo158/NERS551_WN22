#!/usr/bin/python3

# This code correspond to the Part C-a of HW4
# It compares the results of the analytical and the numerical soltion
# run solvePKE and directPKE before running this

#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt


# output file name
input_file_analytical = "PKE_direct.out"
input_file_numerical = "PKE_sol.out"
input_file_sixgroup = "PKE_sol_mg.out"

# get the data from the text file
p_dir = np.loadtxt(input_file_analytical)[:, 1]
p_num = np.loadtxt(input_file_numerical)[:, 1]
p_six = np.loadtxt(input_file_sixgroup)[:, 1]
t = np.loadtxt(input_file_numerical)[:, 3]

# show the two figures
fig, axs = plt.subplots(1, 1)
axs.plot(t, p_dir, t, p_num, t, p_six)
axs.set_xlim(0, 6)
axs.set_xlabel('time')
axs.legend([ 'Analytical','Numerical(1-group)','Numerical (6-group)'])
axs.grid(True)
# plt.show()
fig.savefig("./figures/PKE_sol_validation.png")

# relative difference -- enable if you want to check the results
# rel_diff =  abs( 100* (p_num - p_dir)/ p_dir) 
# axs[1].plot(t,rel_diff)
# axs[1].grid(True)
# axs[1].set_ylabel('crelative difference (%)')
# fig.tight_layout()
# plt.show()
# fig.savefig("./figures/PKE_sol_relDiff.png")




