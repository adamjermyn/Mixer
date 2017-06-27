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

plt.figure(figsize=(6,4.3))

plt.subplot(122)

x, y = delta, r[...,4,3,0]
plt.plot(x,y, label='$v_r v_\\theta$')
x, y = delta, r[...,1,3,0]
plt.plot(x,y, label='$v_r r_\\theta$')
plt.xlabel('$\delta$')
plt.legend(loc='upper left')
plt.tight_layout()

plt.subplot(121)

x, y = delta, r[...,3,3,0]
plt.plot(x,y, label='$v_r v_r$')
x, y = delta, r[...,0,3,0]
plt.plot(x,y, label='$v_r r_r$')
print(flux(np.pi/2,1e-5,'14','Sr'))
plt.plot(x, (1./3)*(flux(np.pi/2, 1e-4, '14', 'Sr')*np.cos(x) + np.sin(x)*flux(np.pi/2, 1e-4, '14', 'St')), label='$v_r r_r$ *')
plt.xlabel('$\delta$')
plt.ylabel('(Mixing Units)')
plt.legend(loc='upper left')
plt.tight_layout()

plt.savefig('Plots/Plot7.pdf')
