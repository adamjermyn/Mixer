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
from multiprocessing import Pool

omega = 10**np.linspace(-3,3,num=100,endpoint=True)

fi = h5py.File('Data/scale_results.dat','w')
fi['omega'] = omega

tS = np.pi/4
tP = np.pi/4
w = 1e-15
tW = np.pi/2
N2 = -1
tolr = 1e-8
tola = 1e-8
maxEval = 10000000

def f(x):
	print(x)
	params = (x, w, tW, tS, tP, N2, tolr, tola, maxEval)
	r = coeffs(params)
	return r

pool = Pool(processes=8)
results = np.array(pool.map(f, omega))

fi['results'] = results
fi.close()
