#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Inputs
y_loc = 10.0
t = 10.0
fname = "iv_profile_t10.csv"
x_0 = -10.0
y_0 = -10.0
g = 1.4
R = 287.0
cv = R/(g-1.0)
alpha = 1.0
rho_inf = 1.0
p_inf = 1.0/g
u_inf = 2.0
v_inf = 2.0
K = 5.0
T_inf = p_inf/(rho_inf*R)
a_inf = np.sqrt(g*R*T_inf)

# Checking to make sure file exists
if not os.path.isfile(fname):
    print("File {} doesn't exist.\nExiting\n\n".format(fname))
    sys.exit()

# Reading data from csv file
data = np.genfromtxt(fname, skip_header=1, delimiter=",")

# Separating data
rho_n = data[:, 0]
rhou_n = data[:, 1]
rhov_n = data[:, 2]
u_n = rhou_n/rho_n
et_n = data[:, 3]
x_n = data[:, 6]

# Constructing exact solution
x_e = np.linspace(-20.0, 20.0, 2000)
xb = x_e - x_0 - u_inf*t
yb = y_loc - y_0 - v_inf*t
rb = np.sqrt(xb**2 + yb**2)
T = T_inf*(1.0 - K**2*(g-1.0)/(8.0*alpha*np.pi**2*a_inf**2)*np.exp(alpha*(1.0-rb**2)))
rho_e = rho_inf*(T/T_inf)**(1.0/(g-1.0))
u_e = u_inf - K/(2.0*np.pi)*yb*np.exp(alpha*(1.0-rb**2)/2.0)
v_e = v_inf + K/(2.0*np.pi)*xb*np.exp(alpha*(1.0-rb**2)/2.0)
et_e = cv*T + 0.5*(u_e**2 + v_e**2)

# Creating comparison plot
plt.plot(x_n, rho_n, "ob", lw=1.4, mfc="None", mec="b")
plt.plot(x_e, rho_e, "-b", lw=1.5)
plt.plot(x_n, et_n, "ok", lw=1.4, mfc="None", mec="k")
plt.plot(x_e, et_e, "-k", lw=1.5)

plt.show()
