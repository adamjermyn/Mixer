#!/usr/bin/env python
import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import sys
sys.path.append(d + '/Python/')

import numpy as np
from pyTurb import correlator
from multiprocessing import Pool

omega = 1.
tS = np.pi/4
tP = np.pi/4
w = 0.01
tW = np.pi/2
N2 = -1
tolr = 1e-30
tola = 1e-30

tRan = np.linspace(0, np.pi, num=200, endpoint=True)
pRan = np.linspace(0, 2*np.pi, num=200, endpoint=True)

correl = np.zeros((len(tRan), len(pRan), 4, 4))
for i,t in enumerate(tRan):
	for j,p in enumerate(pRan):
		correl[i,j] = correlator(1, t, p, 0, 0, 0, omega, w, tW, tS, tP, N2)

correl[np.abs(correl) < 1e-10] = 0

print(correl)

import matplotlib.pyplot as plt
plt.imshow(np.log10(np.abs(correl[...,2,2])))
plt.colorbar()
plt.show()