#!/usr/bin/python3

# This code correspond to the Part C-2 of HW4

import numpy as np

# init parameters
t_max = 6
dt = 1.0E-3
dt_prev = 1.0E-3
n = t_max/dt

LAMBD = 2.6E-15
LAMBD_prev = 2.6E-15
LAMBD0 = 2.6E-15

# one-group beta and lambd
# beta_eff = 0.0076
# beta_eff_prev = 0.0076
# lambd = 0.49405
# lambd_prev = 0.49405

# six-group beta and lambd
beta  = [0.02584,0.152,0.13908, 0.30704, 0.1102, 0.02584]
beta_eff_prev = beta
lambd = [0.0128, 0.0318, 0.119, 0.3181, 1.4027, 3.9286]
lambd_prev = lambd 

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

def k1(x):
    return 0.5-x/6.0+x*x/24.0-x*x*x/120.0+x*x*x*x/720.0

def k0(x):
    return 1.0-x/2.0+x*x/6.0-x*x*x/24.0+x*x*x*x/120.0-x*x*x*x*x/720.0

# output_file = "PKE_sol" + "_dt_" +str("{:.1e}".format(dt)) + "_Lambda_" +str("{:.1e}".format(LAMBD)) + ".out" 
output_file = "PKE_sol.out"
fout = open(output_file,'w')

#w = 1.0
for i in range(int(n)):
    # step 0 -- set the reactivity curve
    t = i/1000
    if i<=1000:
        rho_im = 0.5*beta_eff*t
    else:
        # rho_im = 0.5*beta_eff-0.5*beta_eff*(i-1000)/5000
        rho_im = 0.5*beta_eff-0.1*beta_eff*(t-1)

    # step 1 -- Determine transformation parameter with eq. (23)
    #        -- Prepare xx and xx with eq. (15),(16) and (17),
    #           the exponential integration functions (10)~(13) are called
    alpha = 1.0/dt_prev*np.log(p_prev/p_pprev)
    lambd_tilda = (lambd+alpha)*dt
    OMEGA = LAMBD0/LAMBD*beta_eff*dt*k1(lambd_tilda)
    G_prev = LAMBD0/LAMBD_prev*beta_eff_prev*p_prev
    zeta_hat = np.exp(-lambd*dt)*zeta_prev + np.exp(alpha*dt)*dt*G_prev*(k0(lambd_tilda)-k1(lambd_tilda))

    # step 2 -- Prepare xx, xx and xx with eq. (18),(19)
    tau = lambd*OMEGA
    Sd_hat = lambd*zeta_hat
    Sd_prev = lambd_prev*zeta_prev

    # step 3 -- Prepare a1,b1 with eq. (32) and (33) 
    lambd_H_tilda = (lambd_H+alpha)*dt
    lambd_H_hat = lambd_H*dt
    rho_d = rho-rho_im
    a1 = ffp*gamma_d*dt*k1(lambd_H_tilda)
    b1 = rho_im+np.exp(-lambd_H_hat)*rho_d_prev-P0*gamma_d*dt*k0(lambd_H_hat)+np.exp(alpha*dt)*gamma_d*dt*ffp_prev*p_prev*(k0(lambd_H_tilda)-k1(lambd_H_tilda))

    # step 4 -- Prepare a,b,c with eq. (36), (37), and(38)
    a = theta*dt*a1/LAMBD
    b = theta*dt*(((b1-beta_eff)/LAMBD-alpha)+tau/LAMBD0)-1
    c = theta*dt/LAMBD0*Sd_hat+np.exp(alpha*dt)*((1-theta)*dt*((((rho_prev-beta_eff)/LAMBD_prev-alpha)*p_prev+Sd_prev/LAMBD0)+p_prev))

    # step 5 -- Evaluate the flux level with eq. (39) or (40)
    if a<0:
        p = (-b-np.sqrt(b*b-4*a*c))/(2*a)
    elif a==0:
        p = -c/b
    
    # step 6 -- Check if exponential transformation is acceptable with eq. (24), 
    #           if not, redo step 3and 4,5 with  α=0
    #           Evaluate power level, precursors and reactivity with eq. (22), (14) and (31)
    delta_p_left = np.abs(p - np.exp(alpha*dt)*p_prev)
    delta_p_right = np.abs(p-p_prev-(p_prev-p_pprev)*dt_prev/dt)
    if delta_p_left > delta_p_right:
        alpha = 0.0
        # step 3
        lambd_H_tilda = (lambd_H+alpha)*dt
        lambd_H_hat = lambd_H*dt
        rho_d = rho-rho_im
        a1 = ffp*gamma_d*dt*k1(lambd_H_tilda)
        b1 = rho_im+np.exp(-lambd_H_hat)*rho_d_prev-P0*gamma_d*dt*k0(lambd_H_hat)+np.exp(alpha*dt)*gamma_d*dt*ffp_prev*p_prev*(k0(lambd_H_tilda)-k1(lambd_H_tilda))

        # step 4
        a = theta*dt*a1/LAMBD
        b = theta*dt*(((b1-beta_eff)/LAMBD-alpha)+tau/LAMBD0)-1
        c = theta*dt/LAMBD0*Sd_hat+np.exp(alpha*dt)*((1-theta)*dt*((((rho_prev-beta_eff)/LAMBD_prev-alpha)*p_prev+Sd_prev/LAMBD0)+p_prev))
        
        # step 5
        if a<0:
            p = (-b-np.sqrt(b*b-4*a*c))/(2*a)
        elif a==0:
            p = -c/b
    
    pn.append(p)
    rho_prev = rho
    rho = a1*p+b1
    rhon.append(rho)
    zeta_prev = zeta
    zeta = p*OMEGA + zeta_hat
    zetan.append(zeta)

    fout.write("%12f  %12f  %12f %12f \n" % (rho, p, zeta, i/1000) )

fout.close()