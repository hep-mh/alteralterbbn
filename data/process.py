# numpy
import numpy as np
# matplotlib
import matplotlib.pyplot as plt

hbar = 6.582119514e-25

data = np.loadtxt("sm_cosmo_file_old.dat")

t    = data[:,0]
T    = data[:,1]
dTdt = data[:,2]
Tnu  = data[:,3]
H    = data[:,5]
nb   = data[:,6]

t    *= 1          # s > s
T    *= 1e3        # GeV > MeV
dTdt *= 1e6*hbar   # GeV/s > MeV^2
Tnu  *= 1e3        # GeV > MeV
H    *= 1e3*hbar   # 1/s > MeV
nb   *= 1e9        # GeV^3 > MeV^3

np.savetxt("sm_cosmo_file.dat", np.column_stack([t, T, dTdt, Tnu, H, nb]))