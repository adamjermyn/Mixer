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

B = 10**np.linspace(-2,2,num=60,endpoint=True)

fi = h5py.File('Data/polar_results.dat','w')
fi['B'] = B
fi['pol'] = np.array([[1,0,0], [0,1,0], [0,0,1]])

tB = 0
pB = 0
tS = np.pi/2
tP = np.pi/2
tW = np.pi/2
omega = 0.1
w = 1e-10
N2 = -1.
eps = 1e-20
tolr = 1e-8
tola = 1e-12
maxEval = 100000

results = np.zeros(list(B.shape) + [3,6,6,2])

def f(x):
	print(x)
	
	r = np.zeros([3,6,6,2])

	tB = np.pi
	pB = 0
	params = (x, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps)
	r[1] = coeffs(params, output=output)

	tB = np.pi/2
	pB = 0
	params = (x, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps)
	r[0] = coeffs(params, output=output)

	tB = np.pi/2
	pB = np.pi/2
	params = (x, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEva, eps)
	r[2] = coeffs(params, output=output)

	return r

pool = Pool(processes=8)
results = np.array(pool.map(f, B))

fi['results'] = results
fi.close()
