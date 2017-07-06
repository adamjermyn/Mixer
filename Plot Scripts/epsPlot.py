import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt
plt.style.use('ggplot')

fi = h5py.File('Data/eps_results.dat','r')
eps = np.array(fi['eps'])
r = np.array(fi['results'])
fi.close()

plt.figure(figsize=(6,4.3))

plt.subplot(121)

x, y = 1/eps**0.5, r[...,4,3,0]
y = np.abs(y)
plt.loglog(x,y, label='$v_r v_\\theta$')
x, y = 1/eps**0.5, r[...,1,3,0]
y = np.abs(y)
plt.loglog(x,y, label='$v_r r_\\theta$')
plt.xlabel('$L (h)$')
plt.ylabel('(Mixing Units)')
plt.legend(loc='upper left')
plt.tight_layout()

plt.subplot(122)

x, y = 1/eps**0.5, r[...,5,3,0]
print(y)
y = np.abs(y)
plt.loglog(x,y, label='$v_r v_\phi$')
x, y = 1/eps**0.5, r[...,2,3,0]
print(y)
y = np.abs(y)
plt.loglog(x,y, label='$v_r r_\phi$')
plt.xlabel('$L (h)$')
plt.legend(loc='upper left')
plt.tight_layout()

plt.savefig('Plots/eps.pdf')
