import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import sys
sys.path.append(d + '/Python/')

import numpy as np
import h5py
from pyTurb import coeffs

fi = h5py.File('Data/baro_results.dat','w')

delta = np.linspace(-np.pi,np.pi,num=30,endpoint=True)

fi['delta'] = delta

omega = 0
w = 1e-15

tS = np.pi/2
tP = np.pi/2
tW = np.pi/2
N2 = -1
tolr = 1e-10
tola = 1e-10

results = np.zeros(list(delta.shape) + [7,7,2])

for i in range(delta.shape[0]):
	print(i)
	params = (omega, w, tW, tS + delta[i], tP, N2, tolr, tola)
	r = coeffs(params)
	results[i] = r

fi['results'] = results
fi.close()
