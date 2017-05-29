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

omega = 10**np.linspace(-3, 0, endpoint=True, num=20)

fi = h5py.File('Data/alpha_results.dat','w')
fi['omega'] = omega

tS = np.pi/4
tP = np.pi/4
tW = np.pi/2
N2 = -1
tolr = 1e-10
tola = 1e-10
maxEval = 1000000

output = np.zeros((6,6))
output[3,4] = 1
output[3,1] = 1
output[3,5] = 1
output[3,2] = 1

def g(x):
	params = (x, 0, tW, tS, tP, N2, tolr, tola, maxEval)
	r0 = coeffs(params, output=output)
	params = (x, 1e-4*x, tW, tS, tP, N2, tolr, tola, maxEval)
	r1 = coeffs(params, output=output)
	return (r1 - r0) / (1e-4)

def f(x):
	print(x)
	data = []
	w = np.linspace(0, 1e-3, endpoint=True, num=4)

	for j in range(w.shape[0]):
		params = (x, w[j], tW, tS, tP, N2, tolr, tola, maxEval)
		r = coeffs(params, output=output)
		print(r)
		data.append(r)
	data = np.array(data)
	sh = data.shape[1:]
	data = np.reshape(data, (len(w), -1))
	ft, resid, rank, sing, rcond = np.polyfit(w, data, 1, full=True)
	ft = ft[0].reshape(sh)
	return ft

#pool = Pool(processes=8)
#results = np.array(pool.map(f, omega))
results = np.array(map(f, omega))

fi['results'] = results
fi.close()
