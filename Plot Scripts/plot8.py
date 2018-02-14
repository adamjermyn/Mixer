import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

from KRflux import flux

import numpy as np
import h5py
import matplotlib.pyplot as plt
plt.style.use('ggplot')

fi = h5py.File('Data/latitude_results.dat','r')
theta = np.array(fi['theta'])

plt.figure(figsize=(5, 12))

r = np.array(fi['results_0.1'])



plt.subplot(321)

x, y = theta, r[...,3,3,0]
plt.plot(x,y, label='$v_r v_r / L_0^2 |N|^2$')
x, y = theta, r[...,0,3,0]
plt.plot(x,y, label='$v_r r_r / L_0^2 |N|$')
plt.plot(x, -(1./3)*flux(x, 0.2, '14', 'Sr'), '--', label='$v_r r_r / L_0^2 |N|$ *')
plt.xlabel('$\\theta$')
plt.ylim([0,0.35])
plt.legend(loc='upper left')
plt.title('$\Omega = 0.2 |N|$')
plt.tight_layout()

plt.subplot(322)

x, y = theta, r[...,4,3,0]
plt.plot(x,y, label='$v_r v_\\theta / L_0^2 |N|^2$')
x, y = theta, r[...,1,3,0]
plt.plot(x,y, label='$v_r r_\\theta / L_0^2 |N|$')
plt.plot(x, -(1./3)*flux(x, 0.2, '24', 'Sr'), '--', label='$v_r r_\\theta / L_0^2 |N|$ *')
plt.xlabel('$\\theta$')
plt.ylim([-0.04,0.05])
plt.legend(loc='upper left')
plt.title('$\Omega = 0.1 |N|$')
plt.tight_layout()

r = np.array(fi['results_1.0'])

plt.subplot(323)

x, y = theta, r[...,3,3,0]
plt.plot(x,y, label='$v_r v_r / L_0^2 |N|^2$')
x, y = theta, r[...,0,3,0]
plt.plot(x,y, label='$v_r r_r / L_0^2 |N|$')
plt.plot(x, -(1./3)*flux(x, 1.0, '14', 'Sr'), '--', label='$v_r r_r / L_0^2 |N|$ *')
plt.xlabel('$\\theta$')
plt.ylim([0,0.25])
plt.legend(loc='upper left')
plt.title('$\Omega = |N|$')
plt.tight_layout()

plt.subplot(324)

x, y = theta, r[...,4,3,0]
plt.plot(x,y, label='$v_r v_\\theta / L_0^2 |N|^2$')
x, y = theta, r[...,1,3,0]
plt.plot(x,y, label='$v_r r_\\theta / L_0^2 |N|$')
plt.plot(x, -(1./3)*flux(x, 1.0, '24', 'Sr'), '--', label='$v_r r_\\theta / L_0^2 |N|$ *')
plt.xlabel('$\\theta$')
plt.ylim([-0.03,0.06])
plt.legend(loc='upper left')
plt.title('$\Omega = |N|$')
plt.tight_layout()

r = np.array(fi['results_10.0'])

plt.subplot(325)

x, y = theta, r[...,3,3,0]
plt.plot(x,y, label='$v_r v_r / L_0^2 |N|^2$')
x, y = theta, r[...,0,3,0]
plt.plot(x,y, label='$v_r r_r / L_0^2 |N|$')
plt.plot(x, -(1./3)*flux(x, 5.0, '14', 'Sr'), '--', label='$v_r r_r / L_0^2 |N|$ *')
plt.xlabel('$\\theta$')
plt.ylim([0,0.1])
plt.legend(loc='upper left')
plt.title('$\Omega = 5 |N|$')
plt.tight_layout()

plt.subplot(326)

x, y = theta, r[...,4,3,0]
plt.plot(x,y, label='$v_r v_\\theta / L_0^2 |N|^2$')
x, y = theta, r[...,1,3,0]
plt.plot(x,y, label='$v_r r_\\theta / L_0^2 |N|$')
plt.plot(x, -(1./3)*flux(x, 1.0, '24', 'Sr'), '--', label='$v_r r_\\theta / L_0^2 |N|$ *')
plt.xlabel('$\\theta$')
plt.ylim([-0.025,0.025])
plt.legend(loc='upper left')
plt.title('$\Omega = 10 |N|$')
plt.tight_layout()


plt.savefig('Plots/Plot8.pdf', bbox_inches='tight')
fi.close()
