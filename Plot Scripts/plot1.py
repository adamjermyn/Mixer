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

plt.figure(figsize=(6,4.3))

plt.subplot(121)
x, y = omega[omega<=2], r[omega<=2][...,3,3,0]
plt.plot(x,y,label='$v_r v_r$')
plt.xlim([0,1])
x, y = omega[omega<=2], r[omega<=2][...,0,3,0]
plt.plot(x,y,label='$v_r r_r$')
# KR comparison
plt.plot(x, x*flux(np.pi/4, x, '14', 'Sr'), label='$v_r r_r$ *')

plt.xlim([0,1])
plt.xlabel('$ \Omega/|N|$')
plt.ylabel('(Mixing Units)')
plt.legend(loc='upper right')

plt.subplot(122)
x, y = omega[(omega>=0.5) & (omega<=1000)], r[(omega>=0.5) & (omega<=1000)][...,3,3,0]
plt.xlim([0,max(omega)])
plt.loglog(x,y,label='$v_r v_r$')
x, y = omega[(omega>=0.5) & (omega<=1000)], r[(omega>=0.5) & (omega<=1000)][...,0,3,0]
plt.xlim([0,max(omega)])
plt.loglog(x,y,label='$v_r r_r$')
plt.xlabel('$ \Omega/|N|$')
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig('Plots/Plot1.pdf')
