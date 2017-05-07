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

omega = 10**np.linspace(-3, 3, endpoint=True, num=50)

fi = h5py.File('Data/alpha_results.dat','w')
fi['omega'] = omega

tS = np.pi/4
tP = np.pi/4
tW = np.pi/2
N2 = -1
tolr = 1e-18
tola = 1e-18
maxEval = 10000000

def f(x):
	print(x)
	data = []
	w = np.linspace(0, 1e-2/x, endpoint=True, num=3)

	for j in range(w.shape[0]):
		params = (x, w[j] * x, tW, tS, tP, N2, tolr, tola, maxEval)
		r = coeffs(params)
		data.append(r)
	data = np.array(data)
	sh = data.shape[1:]
	data = np.reshape(data, (len(w), -1))
	ft, resid, rank, sing, rcond = np.polyfit(w * x, data, 1, full=True)
	ft = ft[0].reshape(sh)
	return ft

pool = Pool(processes=8)
results = np.array(pool.map(f, omega))

fi['results'] = results
fi.close()
