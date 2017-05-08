#!/usr/bin/env python
import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import sys
sys.path.append(d + '/Python/')

import numpy as np
from pyTurb import coeffs
from multiprocessing import Pool

omega = 839.557862

tS = np.pi/4
tP = np.pi/4
w = 0
tW = np.pi/2
N2 = -1
tolr = 1e-30
tola = 1e-30
maxEval = 10000


output = np.zeros((6,6))
output[3,3] = 1
output[0,3] = 1
output[3,4] = 1
output[3,1] = 1
output[3,5] = 1
output[3,2] = 1

params = (omega, w, tW, tS, tP, N2, tolr, tola, maxEval)
r = coeffs(params, output=output)
print(r)

omega = 1000.

params = (omega, w, tW, tS, tP, N2, tolr, tola, maxEval)
r = coeffs(params, output=output)
print(r)
