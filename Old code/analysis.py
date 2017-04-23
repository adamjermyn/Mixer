import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter as gf

def ang(mat):
	m = mat[:,:,1:4,-3:]
	return np.array([m[:,:,0,1] - m[:,:,1,0],m[:,:,0,2] - m[:,:,2,0],m[:,:,2,1] - m[:,:,1,2]])

def visc(mat):
	m = mat[:,:,1:4,:3]
	return (np.sum(np.sum(mat**2,axis=-1),axis=-1)**0.5)/3

x = np.loadtxt('output',delimiter=',')

tw = x[:8*60*50]
bw = x[8*60*50:]

tw = np.reshape(tw,(60,50,8,7))
bw = np.reshape(bw,(101,41,8,7))


vtw = visc(tw)
tw = ang(tw)
bw = ang(bw)

tw = np.sum(tw**2,axis=0)**0.5
bw = np.sum(bw**2,axis=0)**0.5

tw = gf(tw,[1,1],mode='nearest')
bw = gf(bw,[1,1],mode='nearest')
vtw = gf(vtw,[1,1],mode='nearest')

tw = np.log10(tw)
bw = np.log10(bw)
vtw = np.log10(vtw)

tw = gf(tw,[1,1],mode='nearest')
bw = gf(bw,[1,1],mode='nearest')
vtw = gf(vtw,[1,1],mode='nearest')

ax = plt.subplot(111)
extent = [0,np.pi,-1,2]
plt.imshow(tw,origin='lower',extent=extent)
aspect = 1.0
ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
cb = plt.colorbar()
cb.set_label('log (Intrinsic Angular Momentum)')
plt.xlabel('$\\theta$')
plt.ylabel('log $\Omega$')
plt.title('Intrinsic Angular Momentum')
plt.savefig('angularTW.eps')
plt.clf()

ax = plt.subplot(111)
extent=[-3.,-1,-3.,2]
plt.imshow(bw,origin='lower',extent=extent)
aspect = 1.0
print extent,extent[1]-extent[0],extent[3]-extent[2],(extent[1]-extent[0])/(extent[3]-extent[2])
print abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect
ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
cb = plt.colorbar()
cb.set_label('log (Intrinsic Angular Momentum)')
plt.xlabel('log $B$')
plt.ylabel('log $\Omega$')
plt.title('Intrinsic Angular Momentum')
plt.savefig('angularBW.eps')
plt.clf()

ax = plt.subplot(111)
extent = [0,np.pi,0,3]
plt.imshow(vtw,origin='lower',extent=extent)
aspect = 1.0
ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
cb = plt.colorbar()
cb.set_label('log (Angle-Averaged Viscosity)')
plt.xlabel('$\\theta$')
plt.ylabel('log $\Omega$')
plt.title('Angle-Averaged Viscosity')
plt.savefig('viscBW.eps')
plt.clf()
