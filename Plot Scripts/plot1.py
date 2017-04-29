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

plt.subplot(121)
x, y = omega[omega<=2], r[omega<=2][...,4,4,0]
plt.plot(x,y,label='$v_r v_r$')
plt.xlim([0,1])
plt.ylim([0,0.25])
x, y = omega[omega<=2], r[omega<=2][...,0,4,0]
plt.plot(x,y,label='$v_r r_r$')
plt.xlim([0,1])
plt.ylim([0,0.25])
plt.xlabel('$ \Omega/|N|$')
plt.ylabel('(Mixing Units)')
plt.legend(loc='upper right')

plt.subplot(122)
x, y = omega[(omega>=0.5) & (omega<=1000)], r[(omega>=0.5) & (omega<=1000)][...,4,4,0]
plt.xlim([0,max(omega)])
plt.loglog(x,y,label='$v_r v_r$')
x, y = omega[(omega>=0.5) & (omega<=1000)], r[(omega>=0.5) & (omega<=1000)][...,0,4,0]
plt.xlim([0,max(omega)])
plt.loglog(x,y,label='$v_r r_r$')
plt.xlabel('$ \Omega/|N|$')
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig('Plots/Plot1.pdf')
