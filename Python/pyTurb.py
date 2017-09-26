import numpy as np
import ctypes

import os
from os.path import dirname, abspath
from tree import tree
d = dirname(dirname(abspath(__file__)))
d = d + '/Build/turb.so'

_turb = ctypes.CDLL(d)
_turb.coeffs2.argtypes = [ctypes.c_double for _ in range(9)] + [ctypes.c_int]
_turb.coeffs2.restype = ctypes.POINTER(ctypes.c_double)
_turb.coeffs3.argtypes = [ctypes.c_double for _ in range(12)] + [ctypes.c_int]
_turb.coeffs3.restype = ctypes.POINTER(ctypes.c_double)

_turb.coeffs2spherical.argtypes = [ctypes.c_double for _ in range(10)] + [ctypes.c_int]
_turb.coeffs2spherical.restype = ctypes.POINTER(ctypes.c_double)
_turb.coeffs3spherical.argtypes = [ctypes.c_double for _ in range(13)] + [ctypes.c_int]
_turb.coeffs3spherical.restype = ctypes.POINTER(ctypes.c_double)

array_type = ctypes.c_int * 36
_turb.coeffs2sphericalSpecific.argtypes = [ctypes.POINTER(ctypes.c_int)] + [ctypes.c_double for _ in range(10)] + [ctypes.c_int]
_turb.coeffs2sphericalSpecific.restype = ctypes.POINTER(ctypes.c_double)
_turb.coeffs3sphericalSpecific.argtypes = [ctypes.POINTER(ctypes.c_int)] + [ctypes.c_double for _ in range(13)] + [ctypes.c_int]
_turb.coeffs3sphericalSpecific.restype = ctypes.POINTER(ctypes.c_double)

array_type = ctypes.c_int * 36
array_type2 = ctypes.c_double * 2
array_type3 = ctypes.c_double * 3
_turb.coeffs2sphericalSpecificBox.argtypes = [ctypes.POINTER(ctypes.c_double) for _ in range(2)] + [ctypes.POINTER(ctypes.c_int)] + [ctypes.c_double for _ in range(9)] + [ctypes.c_int, ctypes.c_double]
_turb.coeffs2sphericalSpecificBox.restype = ctypes.POINTER(ctypes.c_double)
_turb.coeffs3sphericalSpecificBox.argtypes = [ctypes.POINTER(ctypes.c_double) for _ in range(2)] + [ctypes.POINTER(ctypes.c_int)] + [ctypes.c_double for _ in range(12)] + [ctypes.c_int, ctypes.c_double]
_turb.coeffs3sphericalSpecificBox.restype = ctypes.POINTER(ctypes.c_double)

_turb.correlator.argtypes = [ctypes.c_double for _ in range(13)]
_turb.correlator.restype = ctypes.POINTER(ctypes.c_double)

def correlator(k, kT, kP, B, tB, pB, omega, w, tW, tS, tP, N2, eps=1e-20):
	global _turb

	res = _turb.correlator(ctypes.c_double(k),ctypes.c_double(kT),ctypes.c_double(kP),\
			ctypes.c_double(B),ctypes.c_double(tB),ctypes.c_double(pB),\
			ctypes.c_double(omega),ctypes.c_double(w),ctypes.c_double(tW),\
			ctypes.c_double(tS),ctypes.c_double(tP),ctypes.c_double(N2), ctypes.c_double(eps))

	ret = np.zeros((4,4))
	for i in range(4):
		for j in range(4):
			ret[i,j] = res[4*i + j]

	return ret

def coeffs2(omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps=1e-10):
	global _turb

	res = _turb.coeffs2(ctypes.c_double(tolr),ctypes.c_double(tola),\
			ctypes.c_double(omega),ctypes.c_double(w),ctypes.c_double(tW),\
			ctypes.c_double(tS),ctypes.c_double(tP),ctypes.c_double(N2),ctypes.c_int(maxEval), ctypes.c_double(eps))

	ret = np.zeros((6,6,2))
	for i in range(6):
		for j in range(6):
			ret[i,j,0] = res[2*(6*i + j)]
			ret[i,j,1] = res[2*(6*i + j) + 1]
	return ret


def coeffs3(B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps=1e-10):
	global _turb

	res = _turb.coeffs3(ctypes.c_double(tolr),ctypes.c_double(tola),\
			ctypes.c_double(B),ctypes.c_double(tB),ctypes.c_double(pB),\
			ctypes.c_double(omega),ctypes.c_double(w),ctypes.c_double(tW),\
			ctypes.c_double(tS),ctypes.c_double(tP),ctypes.c_double(N2),\
			ctypes.c_int(maxEval), ctypes.c_double(eps))

	ret = np.zeros((6,6,2))
	for i in range(6):
		for j in range(6):
			ret[i,j,0] = res[2*(6*i + j)]
			ret[i,j,1] = res[2*(6*i + j) + 1]

	return ret

def coeffs2spherical(theta, omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps=1e-10):
	global _turb

	res = _turb.coeffs2spherical(ctypes.c_double(theta), ctypes.c_double(tolr),\
			ctypes.c_double(tola),ctypes.c_double(omega),ctypes.c_double(w),ctypes.c_double(tW),\
			ctypes.c_double(tS),ctypes.c_double(tP),ctypes.c_double(N2),ctypes.c_int(maxEval), ctypes.c_double(eps))

	ret = np.zeros((6,6,2))
	for i in range(6):
		for j in range(6):
			ret[i,j,0] = res[2*(6*i + j)]
			ret[i,j,1] = res[2*(6*i + j) + 1]
	return ret

def coeffs3spherical(theta, B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps=1e-10):
	global _turb

	res = _turb.coeffs3spherical(ctypes.c_double(theta), ctypes.c_double(tolr),\
			ctypes.c_double(tola),ctypes.c_double(B),ctypes.c_double(tB),\
			ctypes.c_double(pB),ctypes.c_double(omega),ctypes.c_double(w),\
			ctypes.c_double(tW),ctypes.c_double(tS),ctypes.c_double(tP),\
			ctypes.c_double(N2),ctypes.c_int(maxEval), ctypes.c_double(eps))

	ret = np.zeros((6,6,2))
	for i in range(6):
		for j in range(6):
			ret[i,j,0] = res[2*(6*i + j)]
			ret[i,j,1] = res[2*(6*i + j) + 1]

	return ret

def coeffs2sphericalSpecific(output, theta, omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps=1e-10):
	global _turb

	res = _turb.coeffs2sphericalSpecific(array_type(*output), ctypes.c_double(theta), ctypes.c_double(tolr),\
			ctypes.c_double(tola),ctypes.c_double(omega),ctypes.c_double(w),ctypes.c_double(tW),\
			ctypes.c_double(tS),ctypes.c_double(tP),ctypes.c_double(N2),ctypes.c_int(maxEval), ctypes.c_double(eps))

	ret = np.zeros((6,6,2))
	for i in range(6):
		for j in range(6):
			ret[i,j,0] = res[2*(6*i + j)]
			ret[i,j,1] = res[2*(6*i + j) + 1]
	return ret

def coeffs3sphericalSpecific(output, theta, B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps=1e-10):
	global _turb

	res = _turb.coeffs3sphericalSpecific(array_type(*output), ctypes.c_double(theta), ctypes.c_double(tolr),\
			ctypes.c_double(tola),ctypes.c_double(B),ctypes.c_double(tB),\
			ctypes.c_double(pB),ctypes.c_double(omega),ctypes.c_double(w),\
			ctypes.c_double(tW),ctypes.c_double(tS),ctypes.c_double(tP),\
			ctypes.c_double(N2),ctypes.c_int(maxEval), ctypes.c_double(eps))

	ret = np.zeros((6,6,2))
	for i in range(6):
		for j in range(6):
			ret[i,j,0] = res[2*(6*i + j)]
			ret[i,j,1] = res[2*(6*i + j) + 1]

	return ret

def coeffs2sphericalSpecificBox(mins, maxs, output, theta, omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps=1e-10):
	global _turb

	res = _turb.coeffs2sphericalSpecificBox(array_type2(*mins), array_type2(*maxs), array_type(*output), ctypes.c_double(theta), ctypes.c_double(tolr),\
			ctypes.c_double(tola),ctypes.c_double(omega),ctypes.c_double(w),ctypes.c_double(tW),\
			ctypes.c_double(tS),ctypes.c_double(tP),ctypes.c_double(N2),ctypes.c_int(maxEval), ctypes.c_double(eps))

	ret = np.zeros((6,6,2))
	for i in range(6):
		for j in range(6):
			ret[i,j,0] = res[2*(6*i + j)]
			ret[i,j,1] = res[2*(6*i + j) + 1]
	return ret

def coeffs3sphericalSpecificBox(mins, maxs, output, theta, B, tB, pB, omega, w, tW, tS, tP, N2, tolr, tola, maxEval, eps=1e-10):
	global _turb

	res = _turb.coeffs3sphericalSpecificBox(array_type3(*mins), array_type3(*maxs), array_type(*output), ctypes.c_double(theta), ctypes.c_double(tolr),\
			ctypes.c_double(tola),ctypes.c_double(B),ctypes.c_double(tB),\
			ctypes.c_double(pB),ctypes.c_double(omega),ctypes.c_double(w),\
			ctypes.c_double(tW),ctypes.c_double(tS),ctypes.c_double(tP),\
			ctypes.c_double(N2),ctypes.c_int(maxEval), ctypes.c_double(eps))

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

		if len(params) == 10:
			params2 = [output2, params[4]] + list(params)
			mins = [0.,0.]
			maxs = [np.pi,2*np.pi]
			co = lambda x: coeffs2sphericalSpecificBox(*x)
		elif len(params) == 13:
			params2 = [output2, params[7]] + list(params)
			mins = [0.,0.,0.]
			maxs = [1.,np.pi,2*np.pi]
			co = lambda x: coeffs3sphericalSpecificBox(*x)
		else:
			raise NotImplementedError('Number of parameters does not match any known specification.')


		def f(x):
			c = None
			if len(params) == 10:
				c = correlator(1., x[0], x[1], 0, 0, 0, params[0], params[1], params[2], params[3], params[4], params[5], params[6])
			elif len(params) == 13:
				c = correlator(x[0], x[1], x[2], params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8], params[9])
			else:
				raise NotImplementedError('Number of parameters does not match any known specification.')


			if np.sum(np.abs(c[2:,2:])) > 1e-13:
				return np.sum(np.abs(c[2:,2:]))
			else:
				return 0


		t = tree(mins, maxs, f)
		t.allSplit(3000)
#		print(len(t.nonzero))

		r = np.zeros((6,6,2))

		# For determining the number of evals
		est = sum([c.mean*c.volume for c in t.nonzero])
		vol = sum([c.volume for c in t.nonzero])

#		print('VOL:',vol/(2*np.pi*np.pi),'EST:',est,'NUM:',len(t.nonzero))

		if vol/(2*np.pi*np.pi) > 0.1:
			params3 = [mins, maxs] + params2
			r = co(params3)
		else:
			for c in t.nonzero:
				params3 = [c.mins, c.maxs] + params2
	#			params3[-1] = 10 + int(params3[-1] * c.volume/vol)
				params3[-2] = 10 + int(params3[-2] * c.mean*c.volume/est)
	#			print(c.mean, c.volume, est, params3[-2])
				res = co(params3)
				r += res

	return r

