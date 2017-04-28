import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt

fi = h5py.File('Data/scale_grad_results.dat','r')
w = np.array(fi['w'])
omega = np.array(fi['omega'])
r = np.array(fi['results'])
fi.close()

plt.figure(figsize=(6,5))

plt.subplot(121)

x, y = w/omega, r[...,4,5,0]
plt.plot(x,y, label='$v_r v_\\theta$')
x, y = w/omega, r[...,4,1,0]
plt.plot(x,y, label='$v_r r_\\theta$')
plt.legend(loc='upper right')
plt.ylim([-0.006,0])
plt.xlabel('$|R\\nabla \ln \Omega |$')
plt.ylabel('(Mixing Units)')
plt.tight_layout()

plt.subplot(122)

x, y = w/omega, r[...,4,6,0]
plt.plot(x,y, label='$v_r v_\\phi$')
x, y = w/omega, r[...,4,2,0]
plt.plot(x,y, label='$v_r r_\\phi$')
plt.legend(loc='upper right')
plt.ylim([-0.03,0])
plt.xlabel('$|R\\nabla \ln \Omega |$')
plt.tight_layout()


plt.savefig('Plots/Plot4.pdf')
