#!/usr/bin/env python

import numpy as np
from gas_dynamics import oblique_relations as obr
from gas_dynamics import normal_relations as nor

theta  = 10.0*np.pi/180.0   # Deflection angle
Minfty = 2.0
Tinfty = 288.0
pinfty = 101.325e3*0.1
g = 1.4
R = 287.0

# Finding freestream quantities
uinfty = Minfty*np.sqrt(g*R*Tinfty)*np.cos(theta)
vinfty = -Minfty*np.sqrt(g*R*Tinfty)*np.sin(theta)
#uinfty = Minfty*np.sqrt(g*R*Tinfty)
#vinfty = 0.0
rhoinfty = pinfty/(R*Tinfty)

print("\nFreestream quantities:")
print("rho_inf = {:3.4f} kg/m^3".format(rhoinfty))
print("u_inf   = {:4.4f} m/s".format(uinfty))
print("v_inf   = {:4.4f} m/s".format(vinfty))
print("p_inf   = {:1.6e} Pa\n".format(pinfty))

# Determing shock angle and post-shock state to set initial conditions
beta = obr.beta(Minfty, theta)
Mn1 = Minfty*np.sin(beta)
Mn2 = nor.mach_post_shock(Mn1)
M2 = Mn2/np.sin(beta-theta)
p2qp1 = obr.p_ratio(Minfty, theta)
T2qT1 = obr.T_ratio(Minfty, theta)
T2 = T2qT1*Tinfty
p2 = p2qp1*pinfty
rho2 = p2/(R*T2)
u2 = np.sqrt(g*R*T2)*M2
v2 = 0.0

print("post-shock state:")
print("theta = {:3.3f} deg.".format(theta*180.0/np.pi))
print("betar = {:3.3f} deg.".format((beta-theta)*180.0/np.pi))
print("rho_2 = {:3.4f} kg/m^3".format(rho2))
print("u_2   = {:4.4f} m/s".format(u2))
print("v_2   = {:4.4f} m/s".format(v2))
print("p_2   = {:1.6e} Pa\n".format(p2))

