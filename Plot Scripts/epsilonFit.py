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


fi = h5py.File('Data/eps_results.dat','w')

tS = np.pi/2
tP = np.pi/2
tW = np.pi/4
N2 = 1.
omega = 1e-4
w = 0.1
tolr = 1e-10
tola = 1e-20
maxEval = 3000000

output = np.zeros((6,6))
output[4,3] = 1
output[1,3] = 1
output[5,3] = 1
output[2,3] = 1

eps = 10**np.linspace(-20, -2, num=40, endpoint=True)
fi['eps'] = eps

def f(x):
	params = (omega, w, tW, tS, tP, N2, tolr, tola, maxEval, x)
	r = coeffs(params, output=output)
	print(x)
	print(r)
	return r

pool = Pool(processes=4)
results = np.array(pool.map(f, eps))
#results = np.array(map(f, eps))

fi['results'] = results
fi.close()
