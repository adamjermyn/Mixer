import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt
plt.style.use('ggplot')

fi = h5py.File('Data/latitude_results.dat','r')
theta = np.array(fi['theta'])

plt.figure(figsize=(6, 14))

r = np.array(fi['results_0.1'])

plt.subplot(321)

x, y = theta, r[...,3,3,0]
plt.plot(x,y, label='$v_r v_r$')
x, y = theta, r[...,0,3,0]
plt.plot(x,y, label='$v_r r_r$')
plt.xlabel('$\\theta$')
plt.ylabel('(Mixing Units)')
plt.legend(loc='upper left')
plt.title('$\Omega = 0.2 |N|$')
plt.tight_layout()

plt.subplot(322)

x, y = theta, r[...,4,3,0]
plt.plot(x,y, label='$v_r v_\\theta$')
x, y = theta, r[...,1,3,0]
plt.plot(x,y, label='$v_r r_\\theta$')
plt.xlabel('$\\theta$')
plt.legend(loc='upper left')
plt.title('$\Omega = 0.1 |N|$')
plt.tight_layout()

r = np.array(fi['results_1.0'])

plt.subplot(323)

x, y = theta, r[...,3,3,0]
plt.plot(x,y, label='$v_r v_r$')
x, y = theta, r[...,0,3,0]
plt.plot(x,y, label='$v_r r_r$')
plt.xlabel('$\\theta$')
plt.ylabel('(Mixing Units)')
plt.legend(loc='upper left')
plt.title('$\Omega = |N|$')
plt.tight_layout()

plt.subplot(324)

x, y = theta, r[...,4,3,0]
plt.plot(x,y, label='$v_r v_\\theta$')
x, y = theta, r[...,1,3,0]
plt.plot(x,y, label='$v_r r_\\theta$')
plt.xlabel('$\\theta$')
plt.legend(loc='upper left')
plt.title('$\Omega = |N|$')
plt.tight_layout()

r = np.array(fi['results_10.0'])

plt.subplot(325)

x, y = theta, r[...,3,3,0]
plt.plot(x,y, label='$v_r v_r$')
x, y = theta, r[...,0,3,0]
plt.plot(x,y, label='$v_r r_r$')
plt.xlabel('$\\theta$')
plt.ylabel('(Mixing Units)')
plt.legend(loc='upper left')
plt.title('$\Omega = 5 |N|$')
plt.tight_layout()

plt.subplot(326)

x, y = theta, r[...,4,3,0]
plt.plot(x,y, label='$v_r v_\\theta$')
x, y = theta, r[...,1,3,0]
plt.plot(x,y, label='$v_r r_\\theta$')
plt.xlabel('$\\theta$')
plt.legend(loc='upper left')
plt.title('$\Omega = 10 |N|$')
plt.tight_layout()


plt.savefig('Plots/Plot8.pdf')
fi.close()
