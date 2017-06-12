import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt
plt.style.use('ggplot')

fi = h5py.File('Data/polar_results.dat','r')
B = np.array(fi['B'])
r = np.array(fi['results'])
fi.close()

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
x, y = B, np.abs(r[...,0,3,3,0])
ax.plot(x,y,label='$r$')
x, y = B, np.abs(r[...,1,3,3,0])
ax.plot(x,y,label='$\\theta$')
x, y = B, np.abs(r[...,2,3,3,0])
ax.plot(x,y,label='$\phi$')
ax.set_xscale('log')
ax.set_yscale('log')
plt.legend()
plt.xlabel('$v_A/l$')
plt.ylabel('$|\langle v_r v_r \\rangle|$ (Mixing units)')
plt.tight_layout()
plt.savefig('Plots/vrvr_B.pdf')
plt.clf()

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
x, y = B, np.abs(r[...,0,3,5,0])
ax.plot(x,y,label='$r$')
x, y = B, np.abs(r[...,1,3,5,0])
ax.plot(x,y,label='$\\theta$')
x, y = B, np.abs(r[...,2,3,5,0])
ax.plot(x,y,label='$\phi$')
ax.set_xscale('log')
ax.set_yscale('log')
plt.legend()
plt.xlabel('$v_A/l$')
plt.ylabel('$|\langle v_r v_r \\rangle |$ (Mixing units)')
plt.tight_layout()
plt.savefig('Plots/vrvp_B.pdf')