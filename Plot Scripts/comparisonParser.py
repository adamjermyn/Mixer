import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import sys
sys.path.append(d + '/Python/')
from pyTurb import coeffs

import numpy as np
import h5py
import matplotlib.pyplot as plt

plt.style.use('ggplot')

w = 0
tW = np.pi/2
N2 = -1
eps = 1e-20
tolr = 1e-4
tola = 1e-4
maxEval = 1000

output = np.zeros((6,6))
output[3,3] = 1
output[4,4] = 1
output[5,5] = 1

def f(latitude, rotation):
	latitude = np.pi/2 - np.pi*latitude/180
	params = (rotation, w, tW, latitude, latitude, N2, tolr, tola, maxEval, eps)
	r = coeffs(params, output=output)
	return r

def findOmega(latitude, coriolis):
	omega = 0.0
	err = 1.0
	while err > 1e-2:
		res = f(latitude, omega)
		vrms = (res[3,3,0] + res[4,4,0] + res[5,5,0])**0.5
		omega1 = coriolis * vrms / 6 # 6 because that's how many scale heights Chan 2001 uses
		err = abs(omega1 - omega)
		omega = omega1
	return omega

def num(s):
		try:
			return int(s)
		except ValueError:
			return float(s)


# Chan 2001 Table 1

# We don't actually use this table except in parsing Table 2.

print('Parsing Chan Table 1...')

ChanTable1 = []
fi = open('Data/chan1.dat','r')
for line in fi:
	line = line.strip()
	line = line.split(' ')
	line = [float(num(l)) for l in line]
	ChanTable1.append(line)

ChanTable1 = np.array(ChanTable1)
np.savetxt('Data/chan1parsed.dat', ChanTable1)


# Chan 2001 Table 2

print('Parsing Chan Table 2...')

counter = 0
ChanTable2 = []
fi = open('Data/chan2.dat','r')
for line in fi:
	print(line)
	line = line.strip()
	line = line.split(' ')
	line = [float(num(l)) for l in line]

	# Need to rescale. They have d=1 (box size) and h=1./6 (they have 5h in the convective region and 1h in the stable region
	# at the boundary). We're only plotting ratios so that's fine.
	# To fix omega note that Co = Omega * d / v_rm = 6 Omega / v_rm (in units of the scale height). We can scale this to our
	# model by substituting our v_rm, so Omega = (v_rm / 6) * Co. Unlike Kapyla they don't use the non-rotating v_rm though,
	# so we need to calculate our v_rm to scale it.

	line[2] *= np.pi # Get back into radians
	line[2] = 90 - (180 / np.pi) * line[2] # Now turn into degrees of latitude
	# Iterate to find Omega
	omega = findOmega(line[2], ChanTable1[counter][5])

	# Rearrange into the same format as Kapyla 2004
	line = [omega, line[2], line[0], line[1], 0, 0, (line[3]**2+line[4]**2+line[5]**2)**0.5, line[3], line[4], line[5], 0, 0]

	ChanTable2.append(line)
	counter += 1

ChanTable2 = np.array(ChanTable2)
# Average inferred rotation rates for cases which nominally have the same rotation
nominalRots = list(set(ChanTable2[:,3]))
for r in nominalRots:
	ChanTable2[ChanTable2[:,3] == r, 0] = np.average(ChanTable2[ChanTable2[:,3] == r, 0])

np.savetxt('Data/chan2parsed.dat', ChanTable2)


# Chan 2001 Table 3

print('Parsing Chan Table 3...')

counter = 0
ChanTable3 = []
fi = open('Data/chan3.dat','r')
for line in fi:
	print(line)
	line = line.strip()
	line = line.split(' ')
	line = [float(num(l)) for l in line]

	# Need to rescale. They have d=1 (box size) and h=1./6 (they have 5h in the convective region and 1h in the stable region
	# at the boundary). We're only plotting ratios so that's fine.
	# To fix omega note that Co = Omega * d / v_rm = 6 Omega / v_rm (in units of the scale height). We can scale this to our
	# model by substituting our v_rm, so Omega = (v_rm / 6) * Co. Unlike Kapyla they don't use the non-rotating v_rm though,
	# so we need to calculate our v_rm to scale it.

	line[2] *= np.pi # Get back into radians
	line[2] = 90 - (180 / np.pi) * line[2] # Now turn into degrees of latitude
	# Iterate to find Omega
	omega = findOmega(line[2], ChanTable1[counter][5])

	# Rearrange into the same format as Kapyla 2004
	line = [line[1], omega, line[2], line[5], line[4], line[3], ChanTable1[counter][3]]

	# Now we need to multiply by the RMS value of the relevant pair of velocities.
	# These are in Table 2.

	line[3] *= ChanTable2[counter][7]*ChanTable2[counter][9]
	line[4] *= ChanTable2[counter][8]*ChanTable2[counter][9]
	line[5] *= ChanTable2[counter][7]*ChanTable2[counter][8]

	# Next we divide by the total RMS velocity.
	line[3] /= line[-1]**2
	line[4] /= line[-1]**2
	line[5] /= line[-1]**2

	# Finally we eliminate the last entry.
	line = line[:-1]

	ChanTable3.append(line)
	counter += 1

ChanTable3 = np.array(ChanTable3)
# Average inferred rotation rates for cases which nominally have the same rotation
nominalRots = list(set(ChanTable3[:,0]))
for r in nominalRots:
	ChanTable3[ChanTable3[:,0] == r, 1] = np.average(ChanTable3[ChanTable3[:,0] == r, 1])
ChanTable3 = ChanTable3[:,1:]

np.savetxt('Data/chan3parsed.dat', ChanTable3)

# Kapyla Table 1

print('Parsing Kapyla Table 1...')

KapylaTable1 = []

fi = open('Data/kaplyaTable1.dat','r')
for line in fi:
	line = line.strip()
	line = line.split(' ')
	if len(line) > 1:
		co = None
		if '01' in line[0]:
			co = 0.1
		elif '05' in line[0]:
			co = 0.5
		elif 'Co1-' in line[0]:
			co = 1.0
		elif 'Co2' in line[0]:
			co = 2.0
		elif 'Co4' in line[0]:
			co = 4.0
		elif 'Co7' in line[0]:
			co = 7.0
		elif 'Co10' in line[0]:
			co = 10.0
		line[0] = str(0.125*co) 	# This 0.125 is actually 0.45/2/sqrt(2) and comes from their scale height (0.45),
									# their definition of the coriolis number,
									# and the non-rotating RMS velocity in our model Sqrt[0.25]/0.9.
		line = [float(num(l)) for l in line]
		KapylaTable1.append(line)

KapylaTable1 = np.array(KapylaTable1)
np.savetxt('Data/kapyla1parsed.dat', KapylaTable1)

# Kapyla Table 3

print('Parsing Kapyla Table 3...')

KapylaTable3 = []
fi = open('Data/kaplyaTable3.dat','r')
for line in fi:
	line = line.strip()
	line = line.split(' ')
	if len(line) > 1:
		co = None
		if 'Co11' in line[0]:
			co = 1.0
		elif 'Co10' in line[0]:
			co = 10.0
		line[0] = str(0.125*co) 	# This 0.125 is actually 0.45/2/sqrt(2) and comes from their scale height (0.45),
									# their definition of the coriolis number,
									# and the non-rotating RMS velocity in our model Sqrt[0.25]/0.9.
		line = [float(num(l)) for l in line]
		KapylaTable3.append(line)

KapylaTable3 = np.array(KapylaTable3)

# Rescale.
# They define vt as (1/3) d u_t. In their units d = 1 so
vt = (1./3) * KapylaTable3[:,-1]
KapylaTable3[:,2:-1] *= vt[:,np.newaxis]
KapylaTable3[:,2:-1] *= 1e-3 # They give these in units of 1e-3

# Now normalise against the total velocity squared
KapylaTable3[:,2:-1] /= KapylaTable3[:,-1,np.newaxis]**2

# Now chop that off
KapylaTable3 = KapylaTable3[:,:-1]

# Subtract sheared results
KapylaTable3[::2,2:] -= KapylaTable3[1::2,2:]
KapylaTable3 = KapylaTable3[::2]

# Finally correct because theta -> - theta and r -> -r in their coordinate system.
KapylaTable3[:,3] *= -1
KapylaTable3[:,2] *= -1

np.savetxt('Data/kapyla3parsed.dat', KapylaTable3)
