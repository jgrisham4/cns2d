#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# Setting file names
data_files = ["results{}d.dat".format(i) for i in range(1,3)] + ["exact.dat"]

# Importing data
r1d = np.genfromtxt(data_files[0])
r2d = np.genfromtxt(data_files[1])
r2dv = np.genfromtxt("results2d-venkat.dat")
rex = np.genfromtxt(data_files[2])

# Plotting density
fig1 = plt.figure()
plt.plot(r1d[:, 0], r1d[:, 1], "sk", label="1d", ms=7)
plt.plot(r2d[:, 0], r2d[:, 1], "ob", label="2d-barth")
plt.plot(r2dv[:, 0], r2dv[:, 1], ">r", label="2d-venkat")
plt.plot(rex[:, 0], rex[:, 1], "-g", label="Exact")
plt.xlabel(r"$x$ [m]")
plt.ylabel(r"density [kg/m$^3$]")
plt.ylim([0.1,1.1])
plt.legend(loc=3)

# Saving plot
fig1.savefig("shock-tube-comparison.pdf")

plt.show()

