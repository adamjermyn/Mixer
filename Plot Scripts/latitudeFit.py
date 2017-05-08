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

theta = np.linspace(0,np.pi,num=100,endpoint=True)

fi = h5py.File('Data/latitude_results.dat','w')
fi['theta'] = theta

w = 1e-15
tW = np.pi/2
N2 = -1
tolr = 1e-8
tola = 1e-8
maxEval = 2000000

output = np.zeros((6,6))
output[3,3] = 1
output[0,3] = 1
output[4,3] = 1
output[1,3] = 1

omega = 0.1

def f(x):
	params = (omega, w, tW, x, x, N2, tolr, tola, maxEval)
	r = coeffs(params, output=output)
	return r

pool = Pool(processes=4)
results = np.array(pool.map(f, theta))
fi['results_0.1'] = results

omega = 1.0

def g(x):
	params = (omega, w, tW, x, x, N2, tolr, tola, maxEval)
	r = coeffs(params, output=output)
	return r

pool = Pool(processes=4)
results = np.array(pool.map(g, theta))
fi['results_1.0'] = results


omega = 10.

def h(x):
	params = (omega, w, tW, x, x, N2, tolr, tola, maxEval)
	r = coeffs(params, output=output)
	return r

pool = Pool(processes=4)
results = np.array(pool.map(h, theta))
fi['results_10.0'] = results
fi.close()
