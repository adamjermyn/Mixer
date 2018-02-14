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

plt.figure(figsize=(5,5))

ax = plt.subplot(121)

x, y = 1/eps**0.5, r[...,4,3,0]
y = np.abs(y)
plt.loglog(x,y, label='$v_r v_\\theta / L_0^2 |N|^2$')
x, y = 1/eps**0.5, r[...,1,3,0]
y = np.abs(y)
plt.loglog(x,y, label='$v_r r_\\theta / L_0^2 |N|$')
plt.xlabel('$d / L_0$')
ax.set_xticks([1e1, 1e4, 1e7, 1e10])
plt.legend(loc='lower right')
plt.tight_layout()

ax = plt.subplot(122)

x, y = 1/eps**0.5, r[...,5,3,0]
y = np.abs(y)
plt.loglog(x,y, label='$v_r v_\phi / L_0^2 |N|^2$')
x, y = 1/eps**0.5, r[...,2,3,0]
y = np.abs(y)
plt.loglog(x,y, label='$v_r r_\phi / L_0^2 |N|$')
plt.xlabel('$d / L_0$')
ax.set_xticks([1e1, 1e4, 1e7, 1e10])
plt.legend(loc='lower right')
plt.tight_layout()

plt.savefig('Plots/eps.pdf', bbox_inches='tight')
