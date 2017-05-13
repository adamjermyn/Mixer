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

B = 10**np.linspace(-3,3,num=200,endpoint=True)

tB = np.pi/4
pB = 0
tS = np.pi/2
tP = 0
tW = np.pi/2
omega = 1.
w = -3./2
N2 = 0
tolr = 1e-10
tola = 1e-10
maxEval = 100000


fi = h5py.File('Data/accretion_results.dat','w')
fi['B'] = B


def f(x):
	print(x)
	params = (x, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval)
	r = coeffs(params)
	return r

pool = Pool(processes=8)
results = np.array(pool.map(f, B))

fi['results'] = results
fi.close()
