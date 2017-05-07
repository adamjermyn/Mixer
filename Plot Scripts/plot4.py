import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt
plt.style.use('ggplot')

fi = h5py.File('Data/scale_grad_results.dat','r')
w = np.array(fi['w'])
omega = np.array(fi['omega'])
r = np.array(fi['results'])
fi.close()

plt.figure(figsize=(6,4.3))

plt.subplot(121)

x, y = w/omega, r[...,3,4,0]
plt.plot(x,y, label='$v_r v_\\theta$')
x, y = w/omega, r[...,3,1,0]
plt.plot(x,y, label='$v_r r_\\theta$')
plt.legend(loc='upper right')
plt.xlabel('$|R\\nabla \ln \Omega |$')
plt.ylabel('(Mixing Units)')
plt.tight_layout()

plt.subplot(122)

x, y = w/omega, r[...,3,5,0]
plt.plot(x,y, label='$v_r v_\\phi$')
x, y = w/omega, r[...,3,2,0]
plt.plot(x,y, label='$v_r r_\\phi$')
plt.legend(loc='upper right')
plt.xlabel('$|R\\nabla \ln \Omega |$')
plt.tight_layout()


plt.savefig('Plots/Plot4.pdf')

