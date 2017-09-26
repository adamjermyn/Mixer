import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt
plt.style.use('ggplot')

fi = h5py.File('Data/accretion_results.dat','r')
B = np.array(fi['B'])
r = np.array(fi['results'])
fi.close()

fig = plt.figure(figsize=(5,10))

ax = fig.add_subplot(4,1,1)

x, y, yerr = B, r[...,3,5,0], r[...,3,5,1]
ax.plot(x,y,label='$v_r v_\phi (L_0^2 |N|^2)$')
#ax.fill_between(x, y - yerr, y + yerr, alpha=0.5)
x, y, yerr = B, r[...,0,2,0], r[...,0,2,1]
y *= (B**2/(4*np.pi))
ax.plot(x,-y,label='$-B_r B_\phi/\mu_0\\rho (L_0^2 |N|^2)$')
#ax.fill_between(x, y - yerr, y + yerr, alpha=0.5)

ax.set_xscale('log')
ax.set_yscale('symlog', linthreshy=1e-3)

plt.xlabel('$\\frac{B}{\Omega L_0 \sqrt{\mu_0\\rho}}}$')

plt.legend(loc='upper left')


ax = fig.add_subplot(4,1,2)

x, y, yerr = B, r[...,3,3,0], r[...,3,3,1]
ax.plot(x,y,label='$v_r v_r (L_0^2 |N|^2)$')
#ax.fill_between(x, y - yerr, y + yerr, alpha=0.5)
x, y, yerr = B, r[...,0,0,0], r[...,0,0,1]
y *= (B**2/(4*np.pi))
ax.plot(x,y,label='$B_r B_r/\mu_0\\rho (L_0^2 |N|^2)$')
#ax.fill_between(x, y - yerr, y + yerr, alpha=0.5)

ax.set_xscale('log')
ax.set_yscale('symlog', linthreshy=1e-3)

plt.xlabel('$\\frac{B}{\Omega L_0 \sqrt{\mu_0\\rho}}}$')

plt.legend(loc='center left')


ax = fig.add_subplot(4,1,3)

x, y, yerr = B, r[...,3,4,0], r[...,3,4,1]
ax.plot(x,y,label='$v_r v_\\theta (L_0^2 |N|^2)$')
#ax.fill_between(x, y - yerr, y + yerr, alpha=0.5)
x, y, yerr = B, r[...,0,1,0], r[...,0,1,1]
y *= (B**2/(4*np.pi))
ax.plot(x,-y,label='$-B_r B_\\theta/\mu_0\\rho (L_0^2 |N|^2)$')
#ax.fill_between(x, y - yerr, y + yerr, alpha=0.5)

ax.set_xscale('log')
ax.set_yscale('symlog', linthreshy=1e-3)

plt.xlabel('$\\frac{B}{\Omega L_0 \sqrt{\mu_0\\rho}}}$')

plt.legend(loc='upper left')


ax = fig.add_subplot(4,1,4)

x, y, yerr = B, r[...,4,5,0], r[...,4,5,1]
ax.plot(x,y,label='$v_\phi v_\\theta (L_0^2 |N|^2)$')
#ax.fill_between(x, y - yerr, y + yerr, alpha=0.5)
x, y, yerr = B, r[...,1,2,0], r[...,1,2,1]
y *= (B**2/(4*np.pi))
ax.plot(x,-y,label='$-B_\phi B_\\theta/\mu_0\\rho (L_0^2 |N|^2)$')

#ax.fill_between(x, y - yerr, y + yerr, alpha=0.5)


ax.set_xscale('log')
ax.set_yscale('symlog', linthreshy=1e-3)

plt.xlabel('$\\frac{B}{\Omega L_0 \sqrt{\mu_0\\rho}}}$')

plt.legend(loc='bottom left')

plt.tight_layout()
plt.savefig('Plots/Plot9.pdf')
