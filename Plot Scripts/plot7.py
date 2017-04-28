import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt

fi = h5py.File('Data/baro_results.dat','r')
delta = np.array(fi['delta'])
r = np.array(fi['results'])
fi.close()

plt.figure(figsize=(6,5))

plt.subplot(122)

x, y = delta, r[...,5,4,0]
plt.plot(x,y, label='$v_r v_\\theta$')
x, y = delta, r[...,1,4,0]
plt.plot(x,y, label='$v_r r_\\theta$')
plt.ylim([-0.04,0.03])
plt.xlabel('$\delta$')
plt.ylabel('(Mixing Units)')
plt.legend(loc='lower right')
plt.tight_layout()

plt.subplot(121)

x, y = delta, r[...,4,4,0]
plt.plot(x,y, label='$v_r v_r$')
x, y = delta, r[...,0,4,0]
plt.plot(x,y, label='$v_r r_r$')
plt.xlabel('$\delta$')
plt.ylim([-0.10,0.25])
plt.legend(loc='lower right')
plt.tight_layout()

plt.savefig('Plots/Plot7.pdf')
