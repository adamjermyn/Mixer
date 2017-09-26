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

tRan = np.linspace(0, np.pi, num=1000, endpoint=True)
pRan = np.linspace(0, 2*np.pi, num=2000, endpoint=True)

correl = np.zeros((len(tRan), len(pRan), 4, 4))
for i,t in enumerate(tRan):
	for j,p in enumerate(pRan):
		correl[i,j] = correlator(1., t, p, 0, 0, 0, omega, w, tW, tS, tP, N2, eps=1e-20)
#		correl[i,j] -= correlator(1., t, p, 0, 0, 0, omega, w, tW, tS, tP, N2, eps=1e-20)

#correl /= (w*omega)

correl = np.sum(np.sum(np.abs(correl[...,2:,2:])**2, axis=-1), axis=-1)**0.5
correl[np.abs(correl) < 1e-20] = 1e-20
correl = np.log10(np.abs(correl))

import matplotlib.pyplot as plt
plt.pcolormesh(np.log10(np.abs(correl)))
plt.axis('off')
plt.show()
