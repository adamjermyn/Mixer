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

omega = 1
w = omega * 10**np.linspace(-5,0,num=100,endpoint=True)

fi = h5py.File('Data/richardson_results.dat','w')
fi['omega'] = omega
fi['w'] = w


tS = np.pi/2
tP = np.pi/2
tW = np.pi/2
N2 = 1
tolr = 1e-10
tola = 1e-10
maxEval = 100000

def f(x):
	params = (omega, x, tW, tS, tP, N2, tolr, tola, maxEval)
	r = coeffs(params)
	return r

pool = Pool(processes=4)
results = np.array(pool.map(f, w))

fi['results'] = results
fi.close()
