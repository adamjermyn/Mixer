#!/usr/bin/env python
import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import sys
sys.path.append(d + '/Python/')

import numpy as np
from pyTurb import correlator

tS = 0
tP = 0
tW = np.pi/2
omega = 1.
w = -3./2
N2 = 1

tRan = np.linspace(0, np.pi, num=50, endpoint=True)
pRan = np.linspace(0, 2*np.pi, num=50, endpoint=True)

correl = np.zeros((len(tRan), len(pRan), 4, 4))
for i,t in enumerate(tRan):
	for j,p in enumerate(pRan):
		correl[i,j] = correlator(1., t, p, 0, 0, 0, omega, 2*omega*w, tW, tS, tP, N2)
		correl[i,j] -= correlator(1., t, p, 0, 0, 0, omega, omega*w, tW, tS, tP, N2)
#		correl[i,j] = correlator(1., t, p, 0, 0, 0, omega, omega * w, tW, tS, tP, N2)

correl /= (omega * w)
correl[np.abs(correl) < 1e-12] = 0

import matplotlib.pyplot as plt
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
#plt.imshow(np.log10(np.abs(correl[...,3,3])), origin='lower', extent=[3.14,3.32,1.54,1.575])

#plt.subplot(212)
#plt.imshow(np.log10(np.abs(correl[...,2,1])), origin='lower', extent=[0,2*np.pi,0,np.pi])
#plt.colorbar()
#x = np.loadtxt('/Users/adamjermyn/Dropbox/Research/MLT/Turbulence/Plot Scripts/pts', delimiter=',')
#plt.scatter(x[:100000,1], x[:100000,0], c='r')

plt.show()
