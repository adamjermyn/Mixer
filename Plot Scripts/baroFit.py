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

fi = h5py.File('Data/baro_results.dat','w')

delta = np.linspace(-np.pi,np.pi,num=60,endpoint=True)

fi['delta'] = delta

omega = 0
w = 1e-15

tS = np.pi/2
tP = np.pi/2
tW = np.pi/2
N2 = -1
eps = 1e-20
tolr = 1e-12
tola = 1e-12
maxEval = 10000

output = np.zeros((6,6))
output[4,3] = 1
output[1,3] = 1
output[3,3] = 1
output[0,3] = 1

def f(x):
	print(x)
	params = (omega, w, tW, tS + x, tP, N2, tolr, tola, maxEval, eps)
	r = coeffs(params, output=output)
	return r

pool = Pool(processes=4)
results = np.array(pool.map(f, delta))

fi['results'] = results
fi.close()
