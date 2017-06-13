import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt
plt.style.use('ggplot')

fi = h5py.File('Data/alpha_results.dat','r')
omega = np.array(fi['omega'])
r = np.array(fi['results'])
fi.close()

plt.figure(figsize=(6,4.3))

plt.subplot(121)

print(r)

x, y = omega, r[...,3,4,0]
y = np.abs(y)
plt.loglog(x,y, label='$v_r v_\\theta$')
x, y = omega, r[...,3,1,0]
y = np.abs(y)
plt.loglog(x,y, label='$v_r r_\\theta$')
plt.legend(loc='upper left')
plt.xlabel('$\Omega/|N|$')
plt.ylabel('(Mixing Units)')
plt.tight_layout()

plt.subplot(122)

x, y = omega, r[...,3,5,0]
y = np.abs(y)
plt.loglog(x,y, label='$v_r v_\\phi$')
x, y = omega, r[...,3,2,0]
y = np.abs(y)
plt.loglog(x,y, label='$v_r r_\\phi$')
plt.legend(loc='upper left')
plt.xlabel('$\Omega/|N|$')
plt.tight_layout()


plt.savefig('Plots/Plot5.pdf')

