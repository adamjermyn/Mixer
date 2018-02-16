#!/usr/bin/env python
import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import sys
sys.path.append(d + '/Python/')

import numpy as np
from pyTurb import correlator

tS = np.pi/4
tP = np.pi/4
tW = np.pi/2
omega = 0.1
w = 1e-3 * omega
N2 = 1

tRan = np.linspace(0, np.pi, num=400, endpoint=True)
pRan = np.linspace(0, 2*np.pi, num=400, endpoint=True)

correl = np.zeros((len(tRan), len(pRan), 4, 4))
for i,t in enumerate(tRan):
	for j,p in enumerate(pRan):
		correl[i,j] = correlator(1., t, p, 0, 0, 0, omega, w, tW, tS, tP, N2, eps=1e-20, order=1)
#		correl[i,j] -= correlator(1., t, p, 0, 0, 0, omega, w, tW, tS, tP, N2, eps=1e-20)

#correl /= (w*omega)

correl = np.sum(np.sum(np.abs(correl[...,2:,2:])**2, axis=-1), axis=-1)**0.5
correl[np.abs(correl) < 1e-16] = 0

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
fig, ax = plt.subplots(figsize=(4,4))
im = ax.imshow(np.abs(correl), origin='lower', extent=[0,2*np.pi,0,np.pi], norm=matplotlib.colors.LogNorm(), aspect=2)
ax.set_xlabel('$\phi(q)$')
ax.set_ylabel('$\\theta(q)$')
cbar = fig.colorbar(im, ax=ax)
cbar.set_label('Growth Rate Squared/$|N|^2$')
plt.savefig('Plots/growth.pdf', bbox_inches='tight')

exit()

plt.subplot(221)
plt.imshow(np.log10(np.abs(correl[...,2,2])), origin='lower', extent=[0,2*np.pi,0,np.pi])
plt.colorbar()
plt.subplot(222)
plt.imshow(np.log10(np.abs(correl[...,2,3])), origin='lower', extent=[0,2*np.pi,0,np.pi])
plt.colorbar()
plt.subplot(223)
plt.imshow(np.log10(np.abs(correl[...,3,2])), origin='lower', extent=[0,2*np.pi,0,np.pi])
plt.colorbar()
plt.subplot(224)
plt.imshow(np.log10(np.abs(correl[...,3,3])), origin='lower', extent=[0,2*np.pi,0,np.pi])
plt.colorbar()
plt.show()

exit()

plt.subplot(221)
plt.imshow(np.log10(np.abs(correl[...,0,0])), origin='lower', extent=[0,2*np.pi,0,np.pi])
plt.colorbar()
plt.subplot(222)
plt.imshow(np.log10(np.abs(correl[...,0,1])), origin='lower', extent=[0,2*np.pi,0,np.pi])
plt.colorbar()
plt.subplot(223)
plt.imshow(np.log10(np.abs(correl[...,1,0])), origin='lower', extent=[0,2*np.pi,0,np.pi])
plt.colorbar()
plt.subplot(224)
plt.imshow(np.log10(np.abs(correl[...,1,1])), origin='lower', extent=[0,2*np.pi,0,np.pi])
plt.colorbar()
plt.show()

#plt.imshow(np.log10(np.abs(correl[...,3,3])), origin='lower', extent=[3.14,3.32,1.54,1.575])

#plt.subplot(212)
#plt.imshow(np.log10(np.abs(correl[...,2,1])), origin='lower', extent=[0,2*np.pi,0,np.pi])
#plt.colorbar()
#x = np.loadtxt('/Users/adamjermyn/Dropbox/Research/MLT/Turbulence/Plot Scripts/pts', delimiter=',')
#plt.scatter(x[:100000,1], x[:100000,0], c='r')


