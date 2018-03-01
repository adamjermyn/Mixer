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


omega = 10**np.linspace(-3,3,num=200,endpoint=True)
w = 0.1

fi = h5py.File('Data/richardson_results_2.dat','w')
fi['omega'] = omega
fi['w'] = w


tS = np.pi/2
tP = np.pi/2
tW = np.pi/2
N2 = 1.
eps = 1e-20
tolr = 1e-7
tola = 1e-25
maxEval = 100000000

output = np.zeros((6,6))
output[3,3] = 1
output[4,4] = 1
output[5,5] = 1
output[5,3] = 1

def f(x):
	params = (x, w*x, tW, tS, tP, N2, tolr, tola, maxEval, eps)
	r = coeffs(params, output=output)
	print(x)
	print(r)
	return r

pool = Pool(processes=4)
results = np.array(pool.map(f, omega))

fi['results'] = results
fi.close()
