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

tB = np.pi/4
pB = 0
tS = np.pi/2
tP = 0
tW = np.pi/2
omega = 1.
w = -3./2
N2 = 0
eps = 1e-20
tolr = 1e-6
tola = 1e-10
maxEval = 3000000

B = 1e-10

output = np.ones((6,6))
output[3,5] = 1
output[3,3] = 1



params = (omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps)
#params = (B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps)
r = coeffs(params, output=output)

print(r)
print(r[...,0])
