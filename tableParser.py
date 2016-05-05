import numpy as np
import sys
import table

np.set_printoptions(precision=3,linewidth=90)

def parseTable(fname):
	fi = open(fname)

	tS = np.linspace(0,np.pi,num=20)
	tP = np.linspace(0,2*np.pi,num=41)
	omega = 10**np.linspace(-3,3,num=30)

	data = []

	lines = []
	for i,line in enumerate(fi):
		if i%10 > 1 and i%10 < 9:
			lines.append(line)
		if i%10 == 8:
			lines = [l[:-1].split(' ') for l in lines]
			lines = [[float(l) for l in line if l!=''] for line in lines]
			data.append(np.array(lines))
			lines = []

	data = np.array(data)
	print data.shape
	data = np.reshape(data,(len(tS),len(tP),len(omega),7,7),order='C')
	return tS,tP,omega,data

tS,tP,omega,t = parseTable('outM')

from scipy.interpolate import RegularGridInterpolator as rgi

# Make interpolator
r = rgi([-np.cos(tS),tP,np.log10(omega)],t,method='linear')

# Generate test points
n = 100
x = -np.random.rand(n)
y = np.random.rand(n)*2*np.pi
z = np.random.rand(n)*6-3

for i in range(n):
	q = r((x[i],y[i],z[i]))
	p = table.lvlFinder(np.arccos(-x[i]),y[i],10**z[i],-1,atol=1e-2,rtol=1e-2)
	err = np.abs(p-q)
	err -= 1e-2
	err[err<0] = 0
	print '-----'
	print p
	print q
	err /= np.abs(q)
	err /= 1e-2
	print err
	print np.max(err)
	print '-----'

