#!/usr/bin/env python

import numpy as np

# Inputs
g = 1.4
R = 287.0  # J/(kg-K)
M = 3.0

# Determining freestream inputs for CNS2d
p = 101.325e1  # Pa
T = 288.15     # K
a = np.sqrt(g*R*T)
u = M*a
v = 0.0
rho = p/(R*T)

# Printing info to terminal
print("u_inf   = {:3.4f} m/s".format(u))
print("v_inf   = {:3.4f} m/s".format(v))
print("p_inf   = {:1.6e} Pa".format(p))
print("rho_inf = {:1.6e} kg/m^3".format(rho))

