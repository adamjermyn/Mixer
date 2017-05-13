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

omega = 1.
w = omega * 10**np.linspace(-4,-1,num=30,endpoint=True)

fi = h5py.File('Data/richardson_results.dat','w')
fi['omega'] = omega
fi['w'] = w


tS = np.pi/4
tP = np.pi/4
tW = np.pi/4
N2 = 1.
tolr = 1e-10
tola = 1e-10
maxEval = 30000

output = np.zeros((6,6))
output[4,3] = 1
output[1,3] = 1
output[5,3] = 1
output[2,3] = 1

def f(x):
	print(x)
	params = (omega, x, tW, tS, tP, N2, tolr, tola, maxEval)
	r = coeffs(params, output=output)
	return r

pool = Pool(processes=4)
results = np.array(pool.map(f, w))

fi['results'] = results
fi.close()
