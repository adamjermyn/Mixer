#!/usr/bin/env python
import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import sys
sys.path.append(d + '/Python/')

import numpy as np
from pyTurb import correlator

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker

tS = np.pi/4
tP = np.pi/4
tW = np.pi/2
omega = 0.1
w = 1e-3 * omega
N2 = 1.

maxOrder = 16

for q in range(maxOrder):
	tRan = np.linspace(0, np.pi, num=100, endpoint=True)
	pRan = np.linspace(0, 2*np.pi, num=100, endpoint=True)

	correl = np.zeros((len(tRan), len(pRan), 4, 4))
	for i,t in enumerate(tRan):
		for j,p in enumerate(pRan):
			correl[i,j] = correlator(1., t, p, 0, 0, 0, omega, w, tW, tS, tP, N2, eps=1e-20, order=q)

	correl = np.sum(np.sum(np.abs(correl[...,2:,2:])**2, axis=-1), axis=-1)**0.5
	correl[np.abs(correl) < 1e-15] = 0

	ax = plt.subplot(maxOrder/4,4,q+1)
	plt.title(str(q))
	plt.imshow(np.abs(correl), origin='lower', extent=[0,2*np.pi,0,np.pi], norm=matplotlib.colors.LogNorm(), aspect=2)
	ax.set_xlabel('$\phi(q)$')
	ax.set_ylabel('$\\theta(q)$')
	plt.colorbar()

plt.tight_layout()
plt.show()


