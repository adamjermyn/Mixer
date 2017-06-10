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

B = 10**np.linspace(-2,2,num=20,endpoint=True)

fi = h5py.File('Data/polar_results.dat','w')
fi['B'] = B
fi['pol'] = np.array([[1,0,0], [0,1,0], [0,0,1]])

tB = 0
pB = 0
tS = np.pi/4
tP = np.pi/4
tW = np.pi/4
omega = 0.1
w = 0
N2 = -1.
eps = 1e-20
tolr = 1e-35
tola = 1e-35
maxEval = 100000000

output = np.zeros((6,6))
output[3,3] = 1
output[3,5] = 1

results = np.zeros(list(B.shape) + [3,6,6,2])

def f(x):
	print(x)
	
	r = np.zeros([3,6,6,2])

	tB = np.pi/4
	pB = 0
	params = (x, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps)
	r[0] = coeffs(params, output=output)
	print(r[0])

	tB = 3*np.pi/4
	pB = 0
	params = (x, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps)
	r[1] = coeffs(params, output=output)
	print(r[1])

	tB = np.pi/2
	pB = np.pi/2
	params = (x, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps)
	r[2] = coeffs(params, output=output)
	print(r[2])

	return r

pool = Pool(processes=4)
results = np.array(pool.map(f, B))

fi['results'] = results
fi.close()
