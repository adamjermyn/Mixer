import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

from KRflux import flux

import numpy as np
import h5py
import matplotlib.pyplot as plt
plt.style.use('ggplot')

def paperOne(w, omega):
	# x is |R grad ln Omega|
	return (6/np.pi)*(0.5*np.sin(np.pi/4)*(-0.687)*w + omega*(-0.589))*np.sin(np.pi/4)

fi = h5py.File('Data/scale_grad_results.dat','r')
w = np.array(fi['w'])
omega = np.array(fi['omega'])
r = np.array(fi['results'])
fi.close()

print(r[...,3,4,0])

plt.figure(figsize=(6,4.3))

plt.subplot(121)

x, y = w/omega, r[...,3,4,0]
plt.plot(x,y, label='$v_r v_\\theta$')
plt.fill_between(x, y - r[...,3,4,1], y + r[...,3,4,1], alpha=0.5)
x, y = w/omega, r[...,3,1,0]
plt.plot(x,y, label='$v_r r_\\theta$')
plt.fill_between(x, y - r[...,3,1,1], y + r[...,3,1,1], alpha=0.5)
plt.legend(loc='upper right')
plt.xlabel('$|R\\nabla \ln \Omega |$')
plt.ylabel('(Mixing Units)')
#plt.yscale('log')
plt.tight_layout()

plt.subplot(122)

x, y = w/omega, r[...,3,5,0]
plt.plot(x,y, label='$v_r v_\\phi$')
plt.fill_between(x, y - r[...,3,5,1], y + r[...,3,5,1], alpha=0.5)
x, y = w/omega, r[...,3,2,0]
plt.plot(x,y, label='$v_r r_\\phi$')
plt.fill_between(x, y - r[...,3,2,1], y + r[...,3,2,1], alpha=0.5)

plt.plot(x, flux(np.pi/4,1.,'13','Omr') * (x/np.sqrt(2)) \
			+ flux(np.pi/4,1.,'13','Omt') * (x/np.sqrt(2)) \
			+ flux(np.pi/4,1.,'13','zero'), label='$v_r v_\phi$ *')

plt.plot(x,paperOne(w, omega)/100, label='$v_r v_\phi$ **')

plt.legend(loc='upper right')
plt.xlabel('$|R\\nabla \ln \Omega |$')
#plt.yscale('log')
plt.tight_layout()


plt.savefig('Plots/Plot4.pdf')

