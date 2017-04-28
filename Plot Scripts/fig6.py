import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt

fi = h5py.File('Data/scale_results.dat','r')
omega = np.array(fi['omega'])
r = np.array(fi['results'])
fi.close()

x, y = omega[omega<5], r[omega<5][...,4,0,0]
plt.plot(x,y)
plt.title('Mid-Latitude $r r$ Viscosity')
plt.xlabel('$ \Omega/|N|$')
plt.ylabel('$\langle v_r r \\rangle$')
plt.tight_layout()
plt.savefig('Plots/6_vrr_omega.pdf')