# numpy
import numpy as np


a1 = np.loadtxt("v2.2_abundance_file.dat")
a2 = np.loadtxt("../io/v2.2/abundance_file.dat")

d = np.abs(a1[1:,:]-a2[1:,:])/a1[1:,:]

print(d, np.max(d))