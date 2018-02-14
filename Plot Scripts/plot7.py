import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

from KRflux import flux

import numpy as np
import h5py
import matplotlib.pyplot as plt
plt.style.use('ggplot')

fi = h5py.File('Data/baro_results.dat','r')
delta = np.array(fi['delta'])
r = np.array(fi['results'])
fi.close()

plt.figure(figsize=(5.5,5))

plt.subplot(122)

x, y = delta, r[...,4,3,0]
plt.plot(x,y, label='$v_r v_\\theta / L_0^2 |N|^2$')
x, y = delta, r[...,1,3,0]
plt.plot(x,y, label='$v_r r_\\theta / L_0^2 |N|$')
plt.plot(x, -(1./3)*(flux(np.pi/2, 1e-4, '24', 'Sr')*np.cos(x) + np.sin(x)*flux(np.pi/2, 1e-4, '24', 'St')), '--', label='$v_r r_\\theta / L_0^2 |N|$ *')
plt.xlabel('$\delta$')
plt.legend(loc='lower right')
plt.tight_layout()

plt.subplot(121)

x, y = delta, r[...,3,3,0]
plt.plot(x,y, label='$v_r v_r / L_0^2 |N|^2$')
x, y = delta, r[...,0,3,0]
plt.plot(x,y, label='$v_r r_r / L_0^2 |N|$')
print(flux(np.pi/2,1e-5,'14','Sr'))
plt.plot(x, -(1./3)*(flux(np.pi/2, 1e-4, '14', 'Sr')*np.cos(x) + np.sin(x)*flux(np.pi/2, 1e-4, '14', 'St')), '--', label='$v_r r_r / L_0^2 |N|$ *')
plt.xlabel('$\delta$')
plt.ylim([-0.5, 0.5])
plt.legend(loc='lower left')
plt.tight_layout()

plt.savefig('Plots/Plot7.pdf', bbox_inches='tight')
