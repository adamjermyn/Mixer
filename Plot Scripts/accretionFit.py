import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import sys
sys.path.append(d + '/Python/')

import numpy as np
import h5py
from pyTurb import coeffs


B = np.linspace(-3,3,num=30,endpoint=True)

fi = h5py.File('Data/accretion_results.dat','w')
fi['B'] = B

tB = 0
pB = 0
tS = np.pi/2
tP = np.pi/2
tW = 3*np.pi/2
omega = 1.
w = 3./2
N2 = 0
chi = 0
tolr = 1e-10
tola = 1e-10

results = np.zeros(list(B.shape) + [7,7,2])

for i in range(B.shape[0]):
	print(i)
	params = (B[i], tB, pB, omega, w, tW, tS, tP, N2, chi, tolr, tola)
	r = coeffs(params)
	results[i] = r

fi['results'] = results
fi.close()
