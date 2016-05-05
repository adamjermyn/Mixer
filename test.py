import numpy as np

def correlator(m):
	vals, vecs = np.linalg.eig(m)

	vals[np.real(vals)<0] = 0
	print vals,vals**2
	ret = np.zeros(m.shape)
	for i in range(len(m)):
		vecs[:,i] /= (0.001+np.abs(vecs[-1,i])**2 + np.abs(vecs[-2,i])**2)**0.5
		print vecs[:,i]
		ret += np.outer(vecs[:,i],np.conjugate(vecs[:,i]))*np.real(vals[i])**2
	return ret

m = np.zeros((5,5))

m[3,3] = 1
m[2,3] = 1.4
m[3,2] = 5.4
m[1,1] = 1.2
m[2,4] = 1.4
m[0,3] = 1.999
m[3,0] = -1
m[1,4] = -1
m[4,1] = -1

print correlator(m)
