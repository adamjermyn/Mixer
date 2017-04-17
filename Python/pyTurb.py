import numpy as np
import ctypes

_turb = ctypes.CDLL('/home/asj42/Documents/Turbulence/Build/turb.so')
_turb.coeffs3.argtypes = [ctypes.c_double for _ in range(12)]
_turb.coeffs3.restype = ctypes.POINTER(ctypes.c_double)
_turb.coeffs2.argtypes = [ctypes.c_double for _ in range(8)]
_turb.coeffs2.restype = ctypes.POINTER(ctypes.c_double)

def coeffs2(omega, w, tW, tS, tP, N2, tolr, tola):
	global _turb

	res = _turb.coeffs2(ctypes.c_double(tolr),ctypes.c_double(tola),\
			ctypes.c_double(omega),ctypes.c_double(w),ctypes.c_double(tW),\
			ctypes.c_double(tS),ctypes.c_double(tP),ctypes.c_double(N2))

	ret = np.zeros((7,7,2))
	for i in range(7):
		for j in range(7):
			ret[i,j,0] = res[2*(7*i + j)]
			ret[i,j,1] = res[2*(7*i + j) + 1]

	return ret

def coeffs3(B, tB, pB, omega, w, tW, tS, tP, N2, chi, tolr, tola):
	global _turb

	res = _turb.coeffs3(ctypes.c_double(tolr),ctypes.c_double(tola),\
			ctypes.c_double(B),ctypes.c_double(tB),ctypes.c_double(pB),\
			ctypes.c_double(omega),ctypes.c_double(w),ctypes.c_double(tW),\
			ctypes.c_double(tS),ctypes.c_double(tP),ctypes.c_double(N2),\
			ctypes.c_double(chi))

	ret = np.zeros((7,7,2))
	for i in range(7):
		for j in range(7):
			ret[i,j,0] = res[2*(7*i + j)]
			ret[i,j,1] = res[2*(7*i + j) + 1]

	return ret

def cylindricalToSphericalTensor(tensor, theta):
	transform = np.zeros(tensor.shape)
	
	transform[...,0,0] = np.sin(theta)
	transform[...,0,1] = 0
	transform[...,0,2] = np.cos(theta)

	transform[...,1,0] = np.cos(theta)
	transform[...,1,1] = 0
	transform[...,1,2] = -np.sin(theta)

	transform[...,2,0] = 0
	transform[...,2,1] = 1.
	transform[...,2,2] = 0

	return np.einsum('...kl,...mn,...ln->...km',transform,transform,tensor)


def coeffs(params, spherical=True):
	'''
	Returns the results the turbulence code. Automatically selects the
	integration dimension based on the number of supplied parameters.

	Accepts as input an optional argument spherical. If this is True (default)
	the returned quantities are in spherical coordinates, otherwise they
	are in cylindrical coordinates. Note that errors are not transformed.
	'''
	if len(params) == 8:
		r = coeffs2(*params)
		theta = params[4]
		if spherical:
			r[:3,:3,0] = cylindricalToSphericalTensor(r[:3,:3,0], theta)
			r[4:,:3,0] = cylindricalToSphericalTensor(r[4:,:3,0], theta)
			r[4:,4:,0] = cylindricalToSphericalTensor(r[4:,4:,0], theta)
			r[:3,4:,0] = cylindricalToSphericalTensor(r[:3,4:,0], theta)
	elif len(params) == 12:
		r = coeffs3(*params)
		theta = params[7]
		if spherical:
			r[:3,:3,0] = cylindricalToSphericalTensor(r[:3,:3,0], theta)
			r[4:,:3,0] = cylindricalToSphericalTensor(r[4:,:3,0], theta)
			r[4:,4:,0] = cylindricalToSphericalTensor(r[4:,4:,0], theta)
			r[:3,4:,0] = cylindricalToSphericalTensor(r[:3,4:,0], theta)
	else:
		raise NotImplementedError('Number of parameters does not match any known specification.')
	return r

