#!/usr/bin/env python

import numpy as np

theta  = 10.0*np.pi/180.0   # Deflection angle
Minfty = 2.0
Tinfty = 288.0
pinfty = 101.325e3*0.1
g = 1.4
R = 287.0

uinfty = Minfty*np.sqrt(g*R*Tinfty)*np.cos(theta)
vinfty = -Minfty*np.sqrt(g*R*Tinfty)*np.sin(theta)
rhoinfty = pinfty/(R*Tinfty)

print("\nFreestream quantities:")
print("rho_inf = {:3.4f} kg/m^3".format(rhoinfty))
print("u_inf   = {:4.4f} m/s".format(uinfty))
print("v_inf   = {:4.4f} m/s".format(vinfty))
print("p_inf   = {:1.6e} Pa\n".format(pinfty))



