import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import sys
sys.path.append(d + '/Python/')

import numpy as np
import h5py
from pyTurb import coeffs

omega = 1
w = omega * 10**np.linspace(-5,0,num=100,endpoint=True)

fi = h5py.File('Data/richardson_results.dat','w')
fi['omega'] = omega
fi['w'] = w


tS = np.pi/2
tP = np.pi/2
tW = np.pi/2
N2 = 1
tolr = 1e-10
tola = 1e-10

results = np.zeros(list(w.shape) + [7,7,2])

for i in range(w.shape[0]):
	print(i)
	params = (omega, w[i], tW, tS, tP, N2, tolr, tola)
	r = coeffs(params)
	results[i] = r

fi['results'] = results
fi.close()
