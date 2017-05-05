import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import sys
sys.path.append(d + '/Python/')

import numpy as np
import h5py
from pyTurb import coeffs

omega = 10**np.linspace(-3, 3, endpoint=True, num=100)

fi = h5py.File('Data/alpha_results.dat','w')
fi['omega'] = omega

tS = np.pi/4
tP = np.pi/4
tW = np.pi/2
N2 = -1
tolr = 1e-8
tola = 1e-8
maxEval = 1000000

results = np.zeros(list(omega.shape) + [6,6,2])

for i in range(omega.shape[0]):
	print(i)
	data = []
	w = np.linspace(1e-10, 3e-3/omega[i], endpoint=True, num=8)
	
	for j in range(w.shape[0]):
		params = (omega[i], w[j] * omega[i], tW, tS, tP, N2, tolr, tola, maxEval)
		r = coeffs(params)
		data.append(r)
	data = np.array(data)
	sh = data.shape[1:]
	data = np.reshape(data, (len(w), -1))
	ft, resid, rank, sing, rcond = np.polyfit(w, data, 1, full=True)
	ft = ft[0].reshape(sh)
	print(resid**2/data**2)
	results[i] = ft

fi['results'] = results
fi.close()
