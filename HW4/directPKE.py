#!/usr/bin/python3

# This code correspond to the Part A of HW4

from ast import Lambda
import numpy as np

# init parameters
t_max = 6
dt = 1.0E-3
dt_prev = 1.0E-3
t1 = 1
n = t_max/dt

LAMBD = 2.6E-15
LAMBD_prev = 2.6E-15
LAMBD0 = 2.6E-15

beta_eff = 0.0076
beta_eff_prev = 0.0076
lambd = 0.49405
lambd_prev = 0.49405

alpha = 0
lambd_H = 0
gamma_d = 0
rho_im = 0  # imposed reactivity
rho = 0
rho_prev = 0
rho_d = 0.0
rho_d_prev = 0.0
ffp = 1.0
ffp_prev = 1.0
P0 = 1.0
p = P0*ffp
p_prev = P0*ffp_prev
p_pprev = P0*ffp_prev
zeta = 1.0
zeta_prev = 1.0
theta = 0.5
delta_p_left = 0.0
delta_p_right = 0.0


pn = []
rhon = []
zetan = []

def h(x):
    return -x-np.log(1-x)

# output_file = "PKE_direct" + "_dt_" +str("{:.1e}".format(dt)) + "_Lambda_" +str("{:.1e}".format(LAMBD)) + ".out" 
output_file = "PKE_direct.out"
fout = open(output_file,'w')

#w = 1.0
for i in range(int(n)):
   if i<=1000:
      t = i/1000
      rho_im = 0.5*beta_eff*t
      rho_drv = 0.5*beta_eff
      a = 0.5*beta_eff
      A = (P0)/(beta_eff-rho_im)
      B = (lambd*beta_eff)/rho_drv * h( (rho_drv*t) /beta_eff )
      p = A*np.exp(B)
   else:
      t = i/1000
      rho_im = 0.5*beta_eff-0.1*beta_eff*(t-1)
      a = 0.5*beta_eff
      rho_drv = -0.1*beta_eff
      ka = ( (0.5+0.1) / 0.5)*t1 - (0.1/0.5)*t
      A = (P0*np.exp(B))/(beta_eff- rho_im )
      B = ( lambd*beta_eff/a) * ( h(a*t1/beta_eff) )
      C = (-lambd*beta_eff/a) * ( h(a*ka/beta_eff) )
      p = A*np.exp(B+C) 

   pn.append(p)

   rho = 1
   rhon.append(rho)

   zeta = 1
   zetan.append(zeta)

   fout.write("%12f  %12f  %12f %12f \n" % (rho, p, zeta, i/1000) )

fout.close()