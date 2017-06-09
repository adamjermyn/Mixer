import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt
plt.style.use('ggplot')

fi = h5py.File('Data/richardson_results.dat','r')
w = np.array(fi['w'])
omega = np.array(fi['omega'])
r = np.array(fi['results'])
fi.close()

plt.figure(figsize=(6,4.3))

plt.subplot(121)

print(w/omega)

x, y = w/omega, r[...,4,3,0]
y = np.abs(y)
print(y)
plt.loglog(x,y, label='$v_r v_\\theta$')
#x, y = w/omega, r[...,1,3,0]
#y = np.abs(y)
#plt.loglog(x,y, label='$v_r r_\\theta$')
plt.xlabel('$|R\\nabla \ln \Omega |$')
plt.ylabel('(Mixing Units)')
plt.legend(loc='upper left')
plt.tight_layout()

plt.subplot(122)

x, y = w/omega, r[...,5,3,0]
y = np.abs(y)
plt.loglog(x,y, label='$v_r v_\phi$')
#x, y = w/omega, r[...,2,3,0]
#y = np.abs(y)
#plt.loglog(x,y, label='$v_r r_\phi$')
plt.xlabel('$|R\\nabla \ln \Omega |$')
plt.legend(loc='upper left')
plt.tight_layout()

plt.savefig('Plots/Plot6.pdf')
