#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

d = np.genfromtxt("residuals.dat")
fig1 = plt.figure()
plt.plot(d[:,0], d[:,1], label="continuity")
plt.plot(d[:,0], d[:,2], label="x-momentum")
plt.plot(d[:,0], d[:,3], label="y-momentum")
plt.plot(d[:,0], d[:,4], label="energy")
plt.legend(loc=1)
plt.grid()
plt.yscale("log")

fig1.savefig("residuals.pdf")

plt.show()

