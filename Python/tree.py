import itertools as it
import numpy as np
#from Queue import PriorityQueue
from queue import PriorityQueue

id_counter = 0

class cube:
	def __init__(self, mins, maxs, func, vals):
		global id_counter
		self.mins = mins
		self.maxs = maxs
		self.func = func
		self.vals = vals
		self.id = id_counter
		id_counter += 1

		self.pts = list(it.product(*[[mins[i], maxs[i]] for i in range(len(mins))]))

		for p in self.pts:
			if p not in vals:
				self.vals[p] = self.func(p)

		self.num = 0
		self.mean = 0
		for p in self.pts:
			self.mean += self.vals[p]
			if self.vals[p] != 0:
				self.num += 1
		self.mean /= len(self.pts)


	@property
	def volume(self):
		return np.product(np.array(self.maxs)-np.array(self.mins))
	
	def split(self, i):
		# splits on dimension i
		newMins1 = np.copy(self.mins)
		newMaxs1 = np.copy(self.maxs)
		newMins2 = np.copy(self.mins)
		newMaxs2 = np.copy(self.maxs)

		newMaxs1[i] = (newMins1[i] + newMaxs1[i])/2
		newMins2[i] = (newMins2[i] + newMaxs2[i])/2

		c1 = cube(newMins1, newMaxs1, self.func, self.vals)
		c2 = cube(newMins2, newMaxs2, self.func, self.vals)

		return c1, c2

	def whichSplit(self):
		# Figures out which dimension to split
		# Just picks the split maximizing the number of
		# nonzero values on one side.
		best = [-1, -1e100]
		for i in range(len(self.mins)):
			n = 0
			m = 0
			for p in self.pts:
				if self.vals[p] != 0:
					if p[i] == self.mins[i]:
						n += 1
					elif p[i] == self.maxs[i]:
						m += 1
			n = max(n, m)
			if n > best[1]:
				best[0] = i
				best[1] = n
			elif n == best[1]:
				# Optimize aspect ratio
				if self.maxs[i] - self.mins[i] > self.maxs[best[0]] - self.mins[best[0]]:
					best[0] = i
					best[1] = n
		return best[0]

class tree:
	def __init__(self, mins, maxs, func):
		# top is a cube

		self.tol = 0
		self.vals = {}

		top = cube(mins, maxs, func, self.vals)

		self.cubes = []
		self.toSplit = PriorityQueue()

		self.cubes.append(top)
		self.toSplit.put((-top.num, -np.product(np.array(maxs)-np.array(mins)), top.id, top))


	@property
	def numEvals(self):
		return len(self.vals)
	
	@property
	def nonzero(self):
		return list([c for c in self.cubes if c.num > 0])	

	@property
	def volume(self):
	    return sum([c.volume for c in self.nonzero])	

	def split(self):
		_,_,_,c = self.toSplit.get()
		i = c.whichSplit()
		c1, c2 = c.split(i)
		if c1.num < 2**len(c.mins) and c1.volume > self.tol or c1.num == 0:
			self.toSplit.put((-c1.num, -c1.volume, c1.id, c1))
		if c2.num < 2**len(c.mins) and c2.volume > self.tol or c2.num == 0:
			self.toSplit.put((-c2.num, -c2.volume, c2.id, c2))
		self.cubes.remove(c)
		self.cubes.append(c1)
		self.cubes.append(c2)

		# The tolerance is set so that we split down to 1% of the total active volume.
		self.tol = 1e-2*self.volume


	def allSplit(self, maxSplits):
		counter = 0
		while not self.toSplit.empty() and counter < maxSplits:
			self.split()
			counter += 1
