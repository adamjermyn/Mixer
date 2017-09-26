import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

from KRflux import flux

import numpy as np
import h5py
import matplotlib.pyplot as plt
plt.style.use('ggplot')

fi = h5py.File('Data/scale_results.dat','r')
omega = np.array(fi['omega'])
r = np.array(fi['results'])
fi.close()

plt.figure(figsize=(4,3.5))

x, y = omega, r[...,3,4,0]
y = np.abs(y)
plt.loglog(x,y,label='$v_r v_\\theta (L_0^2 |N|^2)$')
x, y = omega, r[...,3,1,0]
y = np.abs(y)
plt.loglog(x,y,label='$v_r r_\\theta (L_0^2 |N|)$')
plt.loglog(x, np.abs(-(1./3)*flux(np.pi/4, x, '24', 'Sr')), '--', label='$v_r r_\\theta (L_0^2 |N|)$ *')
plt.xlabel('$ \Omega/|N|$')
plt.legend(loc='bottom right')

plt.tight_layout()
plt.savefig('Plots/Plot2.pdf')
