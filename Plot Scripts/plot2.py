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

plt.figure(figsize=(6,5))

x, y = omega, r[...,3,4,0]
y = np.abs(y)
plt.loglog(x,y,label='$v_r v_\\theta$')
x, y = omega, r[...,3,1,0]
y = np.abs(y)
plt.loglog(x,y,label='$v_r r_\\theta$')
plt.xlabel('$ \Omega/|N|$')
plt.ylabel('(Mixing Units)')
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig('Plots/Plot2.pdf')
