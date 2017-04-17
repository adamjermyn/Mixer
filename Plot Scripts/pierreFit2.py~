import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
from pyTurb import coeffs

omega = 10**np.linspace(-5,2,num=25,endpoint=True)
theta =	np.linspace(0,np.pi,num=20,endpoint=True)
wOverOmega = 10**np.linspace(-6,0,num=20,endpoint=True)

fi = h5py.File('Data/results.dat','w')
fi['omega'] = omega
fi['theta'] = theta
fi['wOverOmega'] = wOverOmega

omega, theta, w = np.meshgrid(omega, theta, wOverOmega, indexing='ij')

w *= omega

tW = np.pi/2
N2 = -1
tolr = 1e-10
tola = 1e-10

results = np.zeros(list(omega.shape) + [7,7,2])

for i in range(omega.shape[0]):
	for j in range(omega.shape[1]):
		for k in range(omega.shape[2]):
			print(i,j,k)
			tS = theta[i,j,k]
			tP = theta[i,j,k]
			params = (omega[i,j,k], w[i,j,k], tW, tS, tP, N2, tolr, tola)
			r = coeffs(params)
			results[i,j,k] = r

fi['results'] = results
fi.close()

