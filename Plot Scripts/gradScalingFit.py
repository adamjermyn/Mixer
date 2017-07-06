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

omega = 1.0
w = omega * np.concatenate((10**np.linspace(-4,-2,num=10,endpoint=False), np.linspace(0.01,1.0 ,num=20,endpoint=True)))

fi = h5py.File('Data/scale_grad_results.dat','w')
fi['omega'] = omega
fi['w'] = w

tS = np.pi/4
tP = np.pi/4
tW = np.pi/2
N2 = -1
eps = 1e-20
tolr = 1e-10
tola = 1e-11
maxEval = 1000000


output = np.zeros((6,6))
output[3,4] = 1
output[3,1] = 1
output[3,5] = 1
output[3,2] = 1

def f(x):
	print(x)
	params = (omega, x, tW, tS, tP, N2, tolr, tola, maxEval, eps)
	r = coeffs(params, output=output)
	print(r)
	return r

pool = Pool(processes=4)
results = np.array(pool.map(f, w))

fi['results'] = results
fi.close()
