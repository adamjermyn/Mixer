#!/usr/bin/env python
import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import sys
sys.path.append(d + '/Python/')

import numpy as np
from pyTurb import correlator

B = 1.
tB = 0
pB = 0
tS = 0
tP = 0
tW = np.pi/2
omega = 1.
w = -3./2
N2 = 1
eps = 1e-20

tRan = np.linspace(0, np.pi, num=400, endpoint=True)
pRan = np.linspace(0, 2*np.pi, num=400, endpoint=True)

correl = np.zeros((len(tRan), len(pRan), 4, 4))
for i,t in enumerate(tRan):
	for j,p in enumerate(pRan):
		correl[i,j] = correlator(1., t, p, B, tB, pB, omega, w, tW, tS, tP, N2, eps=1e-2, order=1)

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
plt.show()
