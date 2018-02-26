import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt
plt.style.use('ggplot')

fi = h5py.File('Data/richardson_results_2.dat','r')
omega = np.array(fi['omega'])
r = np.array(fi['results'])
fi.close()


plt.figure(figsize=(5,5))


x, y = omega, r[...,3,3,0]
y = np.abs(y)
print(y)
plt.loglog(x,y, label='$v_r v_r / L_0^2 |N|^2$')

x, y = omega, r[...,4,4,0]
y = np.abs(y)
print(y)
plt.loglog(x,y, label='$v_\\theta v_\\theta / L_0^2 |N|^2$')

x, y = omega, r[...,5,5,0]
y = np.abs(y)
print(y)
plt.loglog(x,y, label='$v_\phi v_\phi / L_0^2 |N|^2$')

x, y = omega, r[...,5,3,0]
y = np.abs(y)
plt.loglog(x,y, label='$v_r v_\phi / L_0^2 |N|^2$')

plt.xlabel('$\Omega/|N|$')
plt.legend(loc='lower right')
plt.tight_layout()

plt.savefig('Plots/Plot6_2.pdf', bbox_inches='tight')
