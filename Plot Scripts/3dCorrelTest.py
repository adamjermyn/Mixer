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
tB = np.pi/4
pB = 0
tS = np.pi/2
tP = 0
tW = np.pi/2
omega = 1.
w = -3./2
N2 = 0
eps = 1e-20

kRan = np.linspace(1e-3, 1-1e-3, num=40, endpoint=True)
tRan = np.linspace(0, np.pi, num=40, endpoint=True)
pRan = np.linspace(0, 2*np.pi, num=40, endpoint=True)

correl = np.zeros((len(kRan), len(tRan), len(pRan), 4))
for l,k in enumerate(kRan):
	for i,t in enumerate(tRan):
		for j,p in enumerate(pRan):
			correl[l,i,j,0] = k/np.max(kRan)
			correl[l,i,j,1] = t/np.max(tRan)
			correl[l,i,j,2] = p/np.max(pRan)
			correl[l,i,j,3] = abs(correlator(k**(-2./3), t, p, B, tB, pB, omega, w, tW, tS, tP, N2))[3,3]

np.savetxt('correl3D.txt',correl.reshape((-1, 4)))


