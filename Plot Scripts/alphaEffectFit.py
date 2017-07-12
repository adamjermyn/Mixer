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

np.set_printoptions(precision=2)

omega = 10**np.linspace(-4.7, -1.5, endpoint=True, num=12) 

fi = h5py.File('Data/alpha_results.dat','w')
fi['omega'] = omega

tS = np.pi/4
tP = np.pi/4
tW = np.pi/2
N2 = -1
eps = 1e-20
tolr = 1e-15
tola = 1e-15
maxEval = 3000000

output = np.zeros((6,6))

output[3,4] = 1
output[3,1] = 1

output[3,5] = 1
output[3,2] = 1


def g(x):
	params = (x, -1e-2*x, tW, tS, tP, N2, tolr, tola, maxEval, eps)
	r0 = coeffs(params, output=output)
	params = (x, 1e-2*x, tW, tS, tP, N2, tolr, tola, maxEval, eps)
	r1 = coeffs(params, output=output)
	print(x)
	print(r0)
	print(r1)
	return (r1 - r0) / (2 * 1e-2*x)

pool = Pool(processes=4)
results = np.array(pool.map(g, omega))
#results = np.array(map(g, omega))

fi['results'] = results
fi.close()

exit()

params = (1e-5, 1e-9, tW, tS, tP, N2, tolr, tola, maxEval, eps)
r0 = coeffs(params, output=output)
print(r0[...,0])
print(r0[...,1])

params = (1e-5, 1e-7, tW, tS, tP, N2, tolr, tola, maxEval, eps)
r0 = coeffs(params, output=output)
print(r0[...,0])
print(r0[...,1])

params = (1e-5, 1e-5, tW, tS, tP, N2, tolr, tola, maxEval, eps)
r0 = coeffs(params, output=output)
print(r0[...,0])
print(r0[...,1])

params = (0, 1e-9, tW, tS, tP, N2, tolr, tola, maxEval, eps)
r0 = coeffs(params, output=output)
print(r0[...,0])
print(r0[...,1])

params = (0, 1e-7, tW, tS, tP, N2, tolr, tola, maxEval, eps)
r0 = coeffs(params, output=output)
print(r0[...,0])
print(r0[...,1])

params = (0, 1e-5, tW, tS, tP, N2, tolr, tola, maxEval, eps)
r0 = coeffs(params, output=output)
print(r0[...,0])
print(r0[...,1])


exit()


