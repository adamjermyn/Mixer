import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pyTurb import correlator
from tree import *

# Funny behaviour occurs between omega=702 and 730. It appears to be a 
# switching-type behaviour, probably due to the floor on unit vectors...

omega = 730
w = 0
tS = np.pi/4
tP = np.pi/4
tW = np.pi/2
N2 = -1



def f(x):
	c = correlator(1, x[0]*np.pi, x[1]*2*np.pi, 0, 0, 0, omega, w, tW, tS, tP, N2)
	print(c)
	if np.sum(np.abs(c[2:,2:])) > 1e-10:
		return 1.0
	else:
		return 0

mins = [0.45,0]
maxs = [0.55,1]

t = tree(mins, maxs, f)
t.allSplit(3000)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111,aspect='equal')
for c in t.cubes:
	if c.num == 2**len(c.mins):
		ax1.add_patch(patches.Rectangle(c.mins, c.maxs[0]-c.mins[0], c.maxs[1]-c.mins[1], fill=True))
	else:
		ax1.add_patch(patches.Rectangle(c.mins, c.maxs[0]-c.mins[0], c.maxs[1]-c.mins[1], fill=False))
plt.show()
exit()


t = tree(mins, maxs, f)
for i in range(10000):
	if t.toSplit.empty():
		break
	t.split()

	if i%10 == 0:
		print(sum([c.num for c in t.cubes]), sum([c.volume for c in t.cubes if c.num > 0]), sum([c.volume for c in t.cubes if c.num == 2**len(c.mins)]))

		fig1 = plt.figure()
		ax1 = fig1.add_subplot(111,aspect='equal')
		for c in t.cubes:
			if c.num == 2**len(c.mins):
				ax1.add_patch(patches.Rectangle(c.mins, c.maxs[0]-c.mins[0], c.maxs[1]-c.mins[1], fill=True))
			else:
				ax1.add_patch(patches.Rectangle(c.mins, c.maxs[0]-c.mins[0], c.maxs[1]-c.mins[1], fill=False))
		plt.show()