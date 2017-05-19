#!/usr/bin/env python
import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import sys
sys.path.append(d + '/Python/')

import numpy as np
from pyTurb import correlator
from multiprocessing import Pool

omega = 1.
tS = np.pi/4
tP = np.pi/4
tW = np.pi/2
N2 = -1
tolr = 1e-30
tola = 1e-30
t = 1.
p = 1.

w = 0.00001
correl = correlator(1, t, p, 0, 0, 0, omega, w, tW, tS, tP, N2)
print(correl)
print('---')
w = 0.0001
correl = correlator(1, t, p, 0, 0, 0, omega, w, tW, tS, tP, N2)
print(correl)
print('---')
w = 0.001
correl = correlator(1, t, p, 0, 0, 0, omega, w, tW, tS, tP, N2)
print(correl)
print('---')
w = 0.01
correl = correlator(1, t, p, 0, 0, 0, omega, w, tW, tS, tP, N2)
print(correl)
print('---')
