#!/usr/bin/env python
import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import sys
sys.path.append(d + '/Python/')

import numpy as np
import h5py
from pyTurb import coeffs

omega = 10**np.linspace(-3,3,num=100,endpoint=True)

fi = h5py.File('Data/scale_results.dat','w')
fi['omega'] = omega

tS = np.pi/4
tP = np.pi/4
w = 1e-15
tW = np.pi/2
N2 = -1
tolr = 1e-9
tola = 1e-9
maxEval = 10000

results = np.zeros((100, 6, 6, 2))

for i,x in enumerate(omega):
	print(i)
	params = (x, w, tW, tS, tP, N2, tolr, tola, maxEval)
	results[i] = coeffs(params)

fi['results'] = results
fi.close()
