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

plt.figure(figsize=(5,5))

plt.subplot(121)
x, y = omega[omega<=2], r[omega<=2][...,3,3,0]
print(y[0])
plt.plot(x,y,label='$v_r v_r / L_0^2 |N|^2 $')
plt.xlim([0,1])
x, y = omega[omega<=2], r[omega<=2][...,0,3,0]
plt.plot(x,y,label='$v_r r_r / L_0^2 |N|$')
print(y[0])
# KR comparison
plt.plot(x, -(1./3)*flux(np.pi/4, x, '14', 'Sr'), '--', label='$v_r r_r / L_0^2 |N|$ *')

plt.xlim([0,1])
plt.xlabel('$ \Omega/|N|$')
plt.legend(loc='lower left')

plt.subplot(122)
x, y = omega[(omega>=0.5) & (omega<=1000)], r[(omega>=0.5) & (omega<=1000)][...,3,3,0]
plt.xlim([1,max(omega)])
plt.loglog(x,y,label='$v_r v_r / L_0^2 |N|^2$')
x, y = omega[(omega>=0.5) & (omega<=1000)], r[(omega>=0.5) & (omega<=1000)][...,0,3,0]
plt.xlim([1,max(omega)])
plt.loglog(x,y,label='$v_r r_r / L_0^2 |N|$')
plt.loglog(x, np.abs(-(1./3)*np.abs(flux(np.pi/4, x, '14', 'Sr'))), '--', label='$v_r r_r / L_0^2 |N|$ *')
plt.xlabel('$ \Omega/|N|$')
plt.legend(loc='lower left')

plt.tight_layout()
plt.savefig('Plots/Plot1.pdf', bbox_inches='tight')
