import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import sys
sys.path.append(d + '/Python/')

import numpy as np
import h5py
from pyTurb import coeffs

omega = 10**np.linspace(-5,2,num=100,endpoint=True)

fi = h5py.File('Data/scale_results.dat','w')
fi['omega'] = omega

tS = np.pi/4
tP = np.pi/4
w = 1e-15
tW = np.pi/2
N2 = -1
tolr = 1e-10
tola = 1e-10

results = np.zeros(list(omega.shape) + [7,7,2])

for i in range(omega.shape[0]):
	print(i)
	params = (omega[i], w, tW, tS, tP, N2, tolr, tola)
	r = coeffs(params)
	results[i] = r

fi['results'] = results
fi.close()
