#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt


# output file name
output_file = "PKE_sol.out"

# time step
dt = 1.0E-3

rho = np.loadtxt(output_file)[:, 0]
p = np.loadtxt(output_file)[:, 1]
zeta = np.loadtxt(output_file)[:, 2]
t = np.loadtxt(output_file)[:, 3]

# plot rho
fig, ax = plt.subplots()
ax.plot(t,rho)
ax.set(xlabel='time (s)', ylabel='reactivity rho',
       title='reactivity rho v.s. time step')
ax.grid()
fig.savefig("./figures/PKE_sol_rho.png")
# plt.show()

# plot p
fig, ax = plt.subplots()
ax.plot(t,p)
ax.set(xlabel='time (s)', ylabel='power p',
       title='power p v.s. time step')
ax.grid()
fig.savefig("./figures/PKE_sol_p.png")
# plt.show()

# plot zeta
fig, ax = plt.subplots()
ax.plot(t,zeta)
ax.set(xlabel='time (s)', ylabel='zeta',
       title='zeta v.s. time step')
ax.grid()
fig.savefig("./figures/PKE_sol_zeta.png")
# plt.show()