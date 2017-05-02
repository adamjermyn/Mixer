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

B = 10**np.linspace(-3,3,num=100,endpoint=True)

tB = 0
pB = 0
tS = np.pi/2
tP = 0
tW = np.pi/2
omega = 1.
w = -3./2
N2 = 0
tolr = 1e-10
tola = 1e-10
maxEval = 100000

B = 0.1

w = 1.
tW = np.pi/2

params = (B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval)
r = coeffs(params)

print(r[3:,3:,0])

w = -1.
tW = 3*np.pi/2

params = (B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval)
r = coeffs(params)

print(r[3:,3:,0])

w = -1.
tW = np.pi/2

params = (B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval)
r = coeffs(params)

print(r[3:,3:,0])

w = 1.
tW = 3*np.pi/2

params = (B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval)
r = coeffs(params)

print(r[3:,3:,0])

w = 10.05
tW = np.pi/2

params = (B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval)
r = coeffs(params)

print(r[3:,3:,0])

w = -10.05
tW = 3*np.pi/2

params = (B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval)
r = coeffs(params)

print(r[3:,3:,0])

w = -10.05
tW = np.pi/2

params = (B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval)
r = coeffs(params)

print(r[3:,3:,0])

w = 10.05
tW = 3*np.pi/2

params = (B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval)
r = coeffs(params)

print(r[3:,3:,0])


exit()
