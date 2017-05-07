import numpy as np
import ctypes

import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
d = d + '/Build/turb.so'

_turb = ctypes.CDLL(d)
_turb.coeffs2.argtypes = [ctypes.c_double for _ in range(8)] + [ctypes.c_int]
_turb.coeffs2.restype = ctypes.POINTER(ctypes.c_double)
_turb.coeffs3.argtypes = [ctypes.c_double for _ in range(11)] + [ctypes.c_int]
_turb.coeffs3.restype = ctypes.POINTER(ctypes.c_double)

_turb.coeffs2spherical.argtypes = [ctypes.c_double for _ in range(9)] + [ctypes.c_int]
_turb.coeffs2spherical.restype = ctypes.POINTER(ctypes.c_double)
_turb.coeffs3spherical.argtypes = [ctypes.c_double for _ in range(12)] + [ctypes.c_int]
_turb.coeffs3spherical.restype = ctypes.POINTER(ctypes.c_double)

array_type = ctypes.c_int * 36
_turb.coeffs2sphericalSpecific.argtypes = [ctypes.POINTER(ctypes.c_int)] + [ctypes.c_double for _ in range(9)] + [ctypes.c_int]
_turb.coeffs2sphericalSpecific.restype = ctypes.POINTER(ctypes.c_double)
_turb.coeffs3sphericalSpecific.argtypes = [ctypes.POINTER(ctypes.c_int)] + [ctypes.c_double for _ in range(12)] + [ctypes.c_int]
_turb.coeffs3sphericalSpecific.restype = ctypes.POINTER(ctypes.c_double)


def coeffs2(omega, w, tW, tS, tP, N2, tolr, tola, maxEval):
	global _turb

	res = _turb.coeffs2(ctypes.c_double(tolr),ctypes.c_double(tola),\
			ctypes.c_double(omega),ctypes.c_double(w),ctypes.c_double(tW),\
			ctypes.c_double(tS),ctypes.c_double(tP),ctypes.c_double(N2),ctypes.c_int(maxEval))

	ret = np.zeros((6,6,2))
	for i in range(6):
		for j in range(6):
			ret[i,j,0] = res[2*(6*i + j)]
			ret[i,j,1] = res[2*(6*i + j) + 1]
	return ret


def coeffs3(B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval):
	global _turb

	res = _turb.coeffs3(ctypes.c_double(tolr),ctypes.c_double(tola),\
			ctypes.c_double(B),ctypes.c_double(tB),ctypes.c_double(pB),\
			ctypes.c_double(omega),ctypes.c_double(w),ctypes.c_double(tW),\
			ctypes.c_double(tS),ctypes.c_double(tP),ctypes.c_double(N2),\
			ctypes.c_int(maxEval))

	ret = np.zeros((6,6,2))
	for i in range(6):
		for j in range(6):
			ret[i,j,0] = res[2*(6*i + j)]
			ret[i,j,1] = res[2*(6*i + j) + 1]

	return ret

def coeffs2spherical(theta, omega, w, tW, tS, tP, N2, tolr, tola, maxEval):
	global _turb

	res = _turb.coeffs2spherical(ctypes.c_double(theta), ctypes.c_double(tolr),\
			ctypes.c_double(tola),ctypes.c_double(omega),ctypes.c_double(w),ctypes.c_double(tW),\
			ctypes.c_double(tS),ctypes.c_double(tP),ctypes.c_double(N2),ctypes.c_int(maxEval))

	ret = np.zeros((6,6,2))
	for i in range(6):
		for j in range(6):
			ret[i,j,0] = res[2*(6*i + j)]
			ret[i,j,1] = res[2*(6*i + j) + 1]
	return ret

def coeffs3spherical(theta, B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval):
	global _turb

	res = _turb.coeffs3spherical(ctypes.c_double(theta), ctypes.c_double(tolr),\
			ctypes.c_double(tola),ctypes.c_double(B),ctypes.c_double(tB),\
			ctypes.c_double(pB),ctypes.c_double(omega),ctypes.c_double(w),\
			ctypes.c_double(tW),ctypes.c_double(tS),ctypes.c_double(tP),\
			ctypes.c_double(N2),ctypes.c_int(maxEval))

	ret = np.zeros((6,6,2))
	for i in range(6):
		for j in range(6):
			ret[i,j,0] = res[2*(6*i + j)]
			ret[i,j,1] = res[2*(6*i + j) + 1]

	return ret

def coeffs2sphericalSpecific(output, theta, omega, w, tW, tS, tP, N2, tolr, tola, maxEval):
	global _turb

	res = _turb.coeffs2sphericalSpecific(array_type(*output), ctypes.c_double(theta), ctypes.c_double(tolr),\
			ctypes.c_double(tola),ctypes.c_double(omega),ctypes.c_double(w),ctypes.c_double(tW),\
			ctypes.c_double(tS),ctypes.c_double(tP),ctypes.c_double(N2),ctypes.c_int(maxEval))

	ret = np.zeros((6,6,2))
	for i in range(6):
		for j in range(6):
			ret[i,j,0] = res[2*(6*i + j)]
			ret[i,j,1] = res[2*(6*i + j) + 1]
	return ret

def coeffs3sphericalSpecific(output, theta, B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval):
	global _turb

	res = _turb.coeffs3sphericalSpecific(array_type(*output), ctypes.c_double(theta), ctypes.c_double(tolr),\
			ctypes.c_double(tola),ctypes.c_double(B),ctypes.c_double(tB),\
			ctypes.c_double(pB),ctypes.c_double(omega),ctypes.c_double(w),\
			ctypes.c_double(tW),ctypes.c_double(tS),ctypes.c_double(tP),\
			ctypes.c_double(N2),ctypes.c_int(maxEval))

	ret = np.zeros((6,6,2))
	for i in range(6):
		for j in range(6):
			ret[i,j,0] = res[2*(6*i + j)]
			ret[i,j,1] = res[2*(6*i + j) + 1]

	return ret

def coeffs(params, output=None):
	'''
	Returns the results the turbulence code. Automatically selects the
	integration dimension based on the number of supplied parameters.

	Accepts as input an optional argument 'output'. If this is None
	then all quantities are computed and output. Otherwise it must be
	an array of integers of shape (6,6) with 1's in the positions where
	outputs are desired and 0's otherwise.
	'''

	if output is None:
		if len(params) == 9:
			params = [params[4]] + list(params)
			r = coeffs2spherical(*params)
		elif len(params) == 12:
			params = [params[7]] + list(params)
			r = coeffs3spherical(*params)
		else:
			raise NotImplementedError('Number of parameters does not match any known specification.')
	else:
		output2 = [0 for _ in range(36)]
		for i in range(6):
			for j in range(6):
				if output[i,j]:
					output2[6*i + j] = 1
		if len(params) == 9:
			params = [output2, params[4]] + list(params)
			r = coeffs2sphericalSpecific(*params)
		elif len(params) == 12:
			params = [output2, params[7]] + list(params)
			r = coeffs3sphericalSpecific(*params)
		else:
			raise NotImplementedError('Number of parameters does not match any known specification.')

	return r

