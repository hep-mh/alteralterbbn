# numpy
import numpy as np
# matplotlib
import matplotlib.pyplot as plt
# sys
from sys import exit


a1 = np.loadtxt("v2.2_abundance_file.dat")
a2 = np.loadtxt("../io/test/abundance_file.dat")

d = np.abs(a1-a2)/a1

print(d)


exit(0)


d1 = np.loadtxt("../io/sm/cosmo_file.dat")
d2 = np.loadtxt("../io/test/cosmo_file.dat")
d3 = np.loadtxt("acropolis_sm_cosmo_file.dat")

i = 4

plt.loglog(d1[:,1], np.abs(d1[:,i]))
plt.loglog(d2[:,1], np.abs(d2[:,i]), ":")
plt.loglog(d3[:,1], np.abs(d3[:,i]), "--")

plt.show()