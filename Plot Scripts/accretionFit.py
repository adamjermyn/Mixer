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

B = 10**np.linspace(-3,1.7,num=40,endpoint=True)

tB = 0
pB = 0
tS = 0
tP = 0
tW = np.pi/2
omega = 1.
w = -3./2
N2 = 1
eps = 1e-20
tolr = 1e-4
tola = 1e-7
maxEval = 30000000


fi = h5py.File('Data/accretion_results.dat','w')
fi['B'] = B

output = np.zeros((6,6))
output[3,5] = 1
output[3,4] = 1
output[3,3] = 1
output[4,5] = 1

def f(x):
	print(x)
	params = (x, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps)
	r = coeffs(params, output=output)
	print(r)
	return r

pool = Pool(processes=2)
results = np.array(pool.map(f, B))

fi['results'] = results
fi.close()
