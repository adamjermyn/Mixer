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

x, y = np.log10(omega), np.log10(np.abs(r[...,5,3,0]))
plt.plot(x,y)
plt.title('Mid-latitude $r \phi$ stress')
plt.xlabel('$\log \Omega/|N|$')
plt.ylabel('$\log|\langle v_r v_\phi\\rangle|$')
plt.tight_layout()
plt.savefig('Plots/3_vrvp_omega.pdf')

