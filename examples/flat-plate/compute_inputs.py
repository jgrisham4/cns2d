#!/usr/bin/env python

import numpy as np
from gas_dynamics import viscosity as v

# Inputs
g    = 1.4
R    = 287.0
L    = 1.0      # m
Tinf = 288.0    # K
ReL  = 35000.0
Minf = 0.3

# Computing primitive variables
uinf = Minf*np.sqrt(g*R*Tinf)
rhoinf = ReL*v.mu(Tinf)/(uinf*L)
vinf = 0.0
pinf = rhoinf*R*Tinf

# Outputting results
print("rho_infty = {:1.6e} kg/m^3".format(rhoinf))
print("u_infty   = {:4.6f} m/s".format(uinf))
print("v_infty   = {:4.6f} m/s".format(vinf))
print("p_infty   = {:6.4f} Pa".format(pinf))
