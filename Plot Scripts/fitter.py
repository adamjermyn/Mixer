import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
from scipy.optimize import curve_fit

def fit(x, y, func, p0, sigma, kind='Relative'):
	popt, _ = curve_fit(func, x, y, p0=p0, sigma=sigma)
	est = func(x, *popt)
	err = np.abs(y - est)
	if kind == 'Relative':
		err = np.abs((y - est)*2/(1e-10 + np.abs(y) + np.abs(est)))
	return popt, np.max(err), np.average(err), est
	
fi = h5py.File('Data/scale_results.dat','r')
omega = np.array(fi['omega'])
r = np.array(fi['results'])
fi.close()

def func1(x, c0, c1, c2):
	return c0*(1 + c1*x**2)/(1 + c2*x**3)

popt, maxerr, avgerr, est = fit(omega, r[...,3,3,0], func1, [1,1,1], 0.1*r[...,3,3,0])
print(popt, maxerr, avgerr)

popt, maxerr, avgerr, est = fit(omega, r[...,0,3,0], func1, [1,1,1], 0.1*r[...,0,3,0])
print(popt, maxerr, avgerr)

def func2(x, c0, c1):
	return c0*(x**2)/(1 + c1*x**3)

popt, maxerr, avgerr, est = fit(omega, r[...,3,4,0], func2, [1,1], 0.1*r[...,3,4,0])
print(popt, maxerr, avgerr)

popt, maxerr, avgerr, est = fit(omega, r[...,3,1,0], func2, [1,1], 0.1*r[...,3,1,0])
print(popt, maxerr, avgerr)

def func5(x, c0, c1):
	return c0*x/(1 + c1*x**3)

popt, maxerr, avgerr, est = fit(omega, r[...,3,5,0], func5, [1,1], 0.1*r[...,3,5,0])
print(popt, maxerr, avgerr)

popt, maxerr, avgerr, est = fit(omega, r[...,3,2,0], func5, [1,1], 0.1*r[...,3,2,0])
print(popt, maxerr, avgerr)


fi = h5py.File('Data/latitude_results.dat','r')
theta = np.array(fi['theta'])
r2 = np.array(fi['results_0.1'])
r1 = np.array(fi['results_1.0'])
r5 = np.array(fi['results_10.0'])
fi.close()

def func3(x, c0, c1, c2):
	return c0 + c1*np.sin(c2*x)

popt, maxerr, avgerr, est = fit(theta, r2[...,3,3,0], func3, [1,1,1], 1e-10 + 0.1*r2[...,3,3,0])
print(popt, maxerr, avgerr)

popt, maxerr, avgerr, est = fit(theta, r2[...,0,3,0], func3, [1,1,1], 1e-10 + 0.1*r2[...,0,3,0])
print(popt, maxerr, avgerr)

popt, maxerr, avgerr, est = fit(theta, r1[...,3,3,0], func3, [1,1,1], 1e-10 + 0.1*r1[...,3,3,0])
print(popt, maxerr, avgerr)

popt, maxerr, avgerr, est = fit(theta, r1[...,0,3,0], func3, [1,1,1], 1e-10 + 0.1*r1[...,0,3,0])
print(popt, maxerr, avgerr)

popt, maxerr, avgerr, est = fit(theta, r5[...,3,3,0], func3, [1,1,1], 1e-10 + 0.1*r5[...,3,3,0])
print(popt, maxerr, avgerr)

popt, maxerr, avgerr, est = fit(theta, r5[...,0,3,0], func3, [1,1,1], 1e-10 + 0.1*r5[...,0,3,0])
print(popt, maxerr, avgerr)

def func4(x, c1):
	return c1*np.sin(2*x)

popt, maxerr, avgerr, est = fit(theta, r2[...,4,3,0], func4, [1], 0.001, kind='Absolute')
print(popt, maxerr, avgerr)

popt, maxerr, avgerr, est = fit(theta, r2[...,1,3,0], func4, [1], 0.001, kind='Absolute')
print(popt, maxerr, avgerr)

popt, maxerr, avgerr, est = fit(theta, r1[...,4,3,0], func4, [1], 0.001, kind='Absolute')
print(popt, maxerr, avgerr)

popt, maxerr, avgerr, est = fit(theta, r1[...,1,3,0], func4, [1], 0.001, kind='Absolute')
print(popt, maxerr, avgerr)

popt, maxerr, avgerr, est = fit(theta, r5[...,4,3,0], func4, [1], 0.001, kind='Absolute')
print(popt, maxerr, avgerr)

popt, maxerr, avgerr, est = fit(theta, r5[...,1,3,0], func4, [1], 0.001, kind='Absolute')
print(popt, maxerr, avgerr)
