import numpy as np
import pickle

p, d, e = pickle.load(open('/data/vault/asj42/dumpP.dat','r'))

omega = np.linspace(-3,3,num=8,endpoint=True)
tw = np.linspace(0,2*np.pi,num=8,endpoint=True)
w = np.linspace(-3,3,num=8,endpoint=True)
tS = np.linspace(0,2*np.pi,num=81,endpoint=True)
tP = np.linspace(0,np.pi,num=41,endpoint=True)

rans = [omega, tw, w, tS, tP]

p[:,0] = np.log10(p[:,0])
p[:,2] = np.log10(p[:,2])

def coords(x, y):
	c = len(y)*(x-min(y))/(max(y)-min(y))
	c = c.astype(int)
	if min(c) == 1:
		c -= 1
	c[c > len(y) - 1] = len(y) - 1
	return c

bins = np.zeros((len(p),5))

for i in range(5):
	bins[:,i] = coords(p[:,i], rans[i])

data = np.zeros((len(omega),len(tw),len(w),len(tS),len(tP),7,7))
errors = np.zeros((len(omega),len(tw),len(w),len(tS),len(tP),7,7))
indicator = np.zeros((len(omega),len(tw),len(w),len(tS),len(tP)))

for i in range(len(bins)):
	data[bins[i][0],bins[i][1],bins[i][2],bins[i][3],bins[i][4]] = d[i]
	errors[bins[i][0],bins[i][1],bins[i][2],bins[i][3],bins[i][4]] = e[i]
	assert indicator[bins[i][0],bins[i][1],bins[i][2],bins[i][3],bins[i][4]] != 1
	indicator[bins[i][0],bins[i][1],bins[i][2],bins[i][3],bins[i][4]] = 1


	if i%5000==0:
		print(i,bins[i],d[i],e[i])

print(np.sum(indicator)*1.0/indicator.size)

pickle.dump([omega, tw, w, tS, tP, data, errors, indicator],open('/data/vault/asj42/dumpG.dat','w+'))
