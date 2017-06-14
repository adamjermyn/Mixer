import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
from scipy.optimize import curve_fit

def fit(x, y, func, p0):
	popt, _ = curve_fit(func, x, y, p0=p0, sigma=0.1*y)
	est = func(x, *popt)
	err = np.abs((y - est)*2/(np.abs(y) + np.abs(est)))
	return popt, np.max(err), np.average(err), est
	
fi = h5py.File('Data/scale_results.dat','r')
omega = np.array(fi['omega'])
r = np.array(fi['results'])
fi.close()

def func1(x, c0, c1, c2):
	return c0*(1 + c1*x**2)/(1 + c2*x**3)

popt, maxerr, avgerr, est = fit(omega, r[...,3,3,0], func1, [1,1,1])
print(popt, maxerr, avgerr)

popt, maxerr, avgerr, est = fit(omega, r[...,0,3,0], func1, [1,1,1])
print(popt, maxerr, avgerr)

def func2(x, c0, c1):
	return c0*(x**2)/(1 + c1*x**3)

popt, maxerr, avgerr, est = fit(omega, r[...,3,4,0], func2, [1,1])
print(popt, maxerr, avgerr)

popt, maxerr, avgerr, est = fit(omega, r[...,3,1,0], func2, [1,1])
print(popt, maxerr, avgerr)

