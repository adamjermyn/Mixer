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

omega = 100

tS = np.pi/4
tP = np.pi/4
w = 1e-15
tW = np.pi/2
N2 = -1
tolr = 1e-18
tola = 1e-18
maxEval = 10000000


params = (omega, w, tW, tS, tP, N2, tolr, tola, maxEval)
r = coeffs(params)
print(r)
