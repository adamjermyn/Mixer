import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt

fi = h5py.File('Data/accretion_results.dat','r')
B = np.array(fi['B'])
r = np.array(fi['results'])
fi.close()

fig = plt.figure(figsize=(5,4))
ax = fig.add_subplot(1,1,1)

x, y = B, r[...,3,5,0]
y = np.abs(y)
ax.plot(x,y,label='$v_r v_\phi$')
x, y = B, r[...,3,2,0]
y = np.abs(y)
ax.plot(x,y,label='$v_r r_\phi$')
ax.set_xscale('log')
plt.xlabel('$\\frac{B}{\Omega l \mu_0 \sqrt{\\rho}}}$')
plt.ylabel('(Mixing Units)')
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig('Plots/Plot9.pdf')
