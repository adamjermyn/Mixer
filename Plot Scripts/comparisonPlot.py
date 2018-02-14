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

# Setup for Mixer

w = 0
tW = np.pi/2
N2 = -1
eps = 1e-20
tolr = 1e-5
tola = 1e-5
maxEval = 10000

output = np.zeros((6,6))
output[3,3] = 1
output[4,4] = 1
output[5,5] = 1
output[4,5] = 1
output[5,3] = 1
output[4,3] = 1

def g(latitude, rotation):
	latitude = np.pi/2 - np.pi*latitude/180
	params = (rotation, w, tW, latitude, latitude, N2, tolr, tola, maxEval, eps)
	r = coeffs(params, output=output)
#	r = np.zeros((6,6,1))
	return r	

ret = g(0, 0)
v_rms0 = (ret[3,3,0] + ret[4,4,0] + ret[5,5,0])**0.5

def f(latitude, coriolis):
	# This method takes as input the coriolis number, defined as
	# Co = Omega * h / v_rms(omega=0).
	# Therefore (using h = 1)
	# Omega * v_rms(omega=0) * Co

	rotation = v_rms0 * coriolis
	return g(latitude, rotation)

# Parse data

ChanTable1 = np.loadtxt('Data/chan1parsed.dat')
ChanTable2 = np.loadtxt('Data/chan2parsed.dat')
ChanTable3 = np.loadtxt('Data/chan3parsed.dat')
KapylaTable1 = np.loadtxt('Data/kapyla1parsed.dat')
KapylaTable3 = np.loadtxt('Data/kapyla3parsed.dat')

# Make plots

# The first plot of interest mixes the two kinds of plots and attempts to pair the nearest possible rotation rates.
# These do not contain our model, just the data.
# This contains the following panels:

# Kapyla - Chan - Type
# 0.022 - 0.000 - diagonal
# 0.113 - 0.130 - diagonal
# 0.225 - 0.260 - diagonal
# 0.450 - 0.524 - diagonal
# 0.900 - 1.048 - diagonal
# 0.225 - 0.260 - off-diagonal
# 0.225 - 0.524 - off-diagonal
# 2.250 - 1.048 - off-diagonal

compareKdiag = [0.022, 0.113, 0.225, 0.450, 0.900]
compareCdiag = [0.000, 0.130, 0.260, 0.524, 1.048]
compareKoffdiag = [0.225, 0.225, 2.250]
compareCoffdiag = [0.260, 0.524, 1.048]

fig, axes = plt.subplots(len(compareKdiag) + len(compareKoffdiag), 1, figsize=(5,12))


latranK = np.linspace(-90,0,num=10,endpoint=True)
latranC = np.linspace(0,90,num=10,endpoint=True)

for i in range(1, len(compareKdiag) + 1):
	print(i)
	dataSelectedK = KapylaTable1[np.abs(KapylaTable1[:,0] - compareKdiag[i-1]) < 0.03]
	dataSelectedC = ChanTable2[np.abs(ChanTable2[:,0] - compareCdiag[i-1]) < 0.03]


	resK = np.array([f(latitude, compareKdiag[i-1]) for latitude in latranK])
	netVK = resK[:,3,3,0] + resK[:,4,4,0] + resK[:,5,5,0]
	resC = np.array([f(latitude, compareCdiag[i-1]) for latitude in latranC])
	netVC = resC[:,3,3,0] + resC[:,4,4,0] + resC[:,5,5,0]

	a, = axes[i-1].plot(latranK, (resK[:,4,4,0]/netVK)**0.5)
	axes[i-1].plot(latranC, (resC[:,4,4,0]/netVC)**0.5, c=a.get_color())
	axes[i-1].scatter(dataSelectedK[:,1], dataSelectedK[:,7]/dataSelectedK[:,6], c=a.get_color(), label='$v_\\theta/v$ (Kapyla)')
	axes[i-1].scatter(dataSelectedC[:,1], dataSelectedC[:,7]/dataSelectedC[:,6], c=a.get_color(), marker='x', label='$v_\\theta/v$ (Chan)')

	a, = axes[i-1].plot(latranK, (resK[:,3,3,0]/netVK)**0.5)
	axes[i-1].plot(latranC, (resC[:,3,3,0]/netVC)**0.5, c=a.get_color())
	axes[i-1].scatter(dataSelectedK[:,1], dataSelectedK[:,9]/dataSelectedK[:,6], c=a.get_color(), label='$v_r/v$ (Kapyla)')
	axes[i-1].scatter(dataSelectedC[:,1], dataSelectedC[:,9]/dataSelectedC[:,6], c=a.get_color(), marker='x', label='$v_r/v$ (Chan)')

	a, = axes[i-1].plot(latranK, (resK[:,5,5,0]/netVK)**0.5)
	axes[i-1].plot(latranC, (resC[:,5,5,0]/netVC)**0.5, c=a.get_color())
	axes[i-1].scatter(dataSelectedK[:,1], dataSelectedK[:,8]/dataSelectedK[:,6], c=a.get_color(), label='$v_\phi/v$ (Kapyla)')
	axes[i-1].scatter(dataSelectedC[:,1], dataSelectedC[:,8]/dataSelectedC[:,6], c=a.get_color(), marker='x', label='$v_\phi/v$ (Chan)')

	axes[i-1].axvline(x=0.5, ymin=0, ymax=1, c='k')

	axes[i-1].text(0.02,0.05,'$\mathrm{Co}$ = '+str(round(compareKdiag[i-1], 3)), transform=axes[i-1].transAxes)
	axes[i-1].text(0.52,0.05,'$\mathrm{Co}$ = '+str(round(compareCdiag[i-1], 3)), transform=axes[i-1].transAxes)
	axes[i-1].grid(False)
	axes[i-1].xaxis.set_visible(False)

	axes[i-1].set_ylim([-0.1,1.1])
	axes[i-1].set_xlim([-100,100])
	axes[i-1].locator_params(axis='y', nbins=4)

for i in range(len(compareKdiag) + 1, len(compareKdiag) + len(compareKoffdiag) + 1):
	dataSelectedK = KapylaTable3[np.abs(KapylaTable3[:,0] - compareKoffdiag[i-1-len(compareKdiag)]) < 0.03]
	dataSelectedC = ChanTable3[np.abs(ChanTable3[:,0] - compareCoffdiag[i-1-len(compareKdiag)]) < 0.03]

	resK = np.array([f(latitude, compareKoffdiag[i-1-len(compareKdiag)]) for latitude in latranK])
	netVK = resK[:,3,3,0] + resK[:,4,4,0] + resK[:,5,5,0]
	resC = np.array([f(latitude, compareCoffdiag[i-1-len(compareKdiag)]) for latitude in latranC])
	netVC = resC[:,3,3,0] + resC[:,4,4,0] + resC[:,5,5,0]

	print(resK[:,4,3,0])
	print(resC[:,4,3,0])
	print(resK[:,5,3,0])
	print(resC[:,5,3,0])
	print(resK[:,4,5,0])
	print(resC[:,4,5,0])
	print('---')

	a, = axes[i-1].plot(latranK, resK[:,4,3,0]/netVK)
	axes[i-1].plot(latranC, resC[:,4,3,0]/netVC, c=a.get_color())
	axes[i-1].scatter(dataSelectedK[:,1], dataSelectedK[:,4], c=a.get_color(), label='$v_r v_\\theta / v^2 $ (Kapyla)')
	axes[i-1].scatter(dataSelectedC[:,1], dataSelectedC[:,4], c=a.get_color(), marker='x', label='$v_r v_\\theta / v^2 $ (Chan)')

	a, = axes[i-1].plot(latranK, resK[:,5,3,0]/netVK)
	axes[i-1].plot(latranC, resC[:,5,3,0]/netVC, c=a.get_color())
	axes[i-1].scatter(dataSelectedK[:,1], dataSelectedK[:,3], c=a.get_color(), label='$v_r v_\phi / v^2 $ (Kapyla)')
	axes[i-1].scatter(dataSelectedC[:,1], dataSelectedC[:,3], c=a.get_color(), marker='x', label='$v_r v_\phi / v^2 $ (Chan)')

	a, = axes[i-1].plot(latranK, resK[:,4,5,0]/netVK)
	axes[i-1].plot(latranC, resC[:,4,5,0]/netVC, c=a.get_color())
	axes[i-1].scatter(dataSelectedK[:,1], dataSelectedK[:,2], c=a.get_color(), label='$v_\\theta v_\phi / v^2 $ (Kapyla)')
	axes[i-1].scatter(dataSelectedC[:,1], dataSelectedC[:,2], c=a.get_color(), marker='x', label='$v_\\theta v_\phi / v^2 $ (Chan)')

	axes[i-1].axvline(x=0.5, ymin=0, ymax=1, c='k')

	axes[i-1].text(0.02,0.05,'$\mathrm{Co}$ = '+str(round(compareKoffdiag[i-1-len(compareKdiag)], 3)), transform=axes[i-1].transAxes)
	axes[i-1].text(0.52,0.05,'$\mathrm{Co}$ = '+str(round(compareCoffdiag[i-1-len(compareKdiag)], 3)), transform=axes[i-1].transAxes)
	axes[i-1].grid(False)

	axes[i-1].xaxis.set_visible(False)

	axes[i-1].set_ylim([-0.23,0.23])
	axes[i-1].set_xlim([-100,100])
	axes[i-1].locator_params(axis='y', nbins=4)

axes[-1].xaxis.set_visible(True)
axes[-1].set_xlabel('Latitude (degrees)')

h, l = axes[-1].get_legend_handles_labels()
# Interleave correctly
h = [h[2*i] for i in range(3)] + [h[2*i+1] for i in range(3)]
l = [l[2*i] for i in range(3)] + [l[2*i+1] for i in range(3)]
legend = fig.legend(h, l, ncol=2, bbox_to_anchor=[0.91, 0.095])

fig.subplots_adjust(wspace=0, hspace=0, bottom=0.13)

plt.savefig('Plots/DataComparison.pdf',bbox_inches='tight')

# The second plot compares <v_i^2>/<v^2> between our model and both simulations.
# This makes use of Kapyla's Table 1 and Chan's Table 2.

data = np.array(list(KapylaTable1) + list(ChanTable2))

rots = list(set(data[:,0]))
rots.sort()


latran = np.linspace(-90, 90, num=15, endpoint=True)

fig, axes = plt.subplots(len(rots)//2, 2, figsize=(10,12))
 
axes = [b for a in axes for b in a]

for i in range(1,len(rots) + 1):
	dataSelected = data[np.abs(data[:,0] - rots[i-1]) < 1e-3]

	res = np.array([f(latitude, rots[i-1]) for latitude in latran])
	netV = res[:,3,3,0] + res[:,4,4,0] + res[:,5,5,0]


	a1, = axes[i-1].plot(latran, (res[:,4,4,0]/netV)**0.5, label='$v_\\theta/v$')
	axes[i-1].scatter(dataSelected[:,1], dataSelected[:,7]/dataSelected[:,6], c=a1.get_color())

	a1, = axes[i-1].plot(latran, (res[:,3,3,0]/netV)**0.5, label='$v_r/v$')
	axes[i-1].scatter(dataSelected[:,1], dataSelected[:,9]/dataSelected[:,6], c=a1.get_color())

	a1, = axes[i-1].plot(latran, (res[:,5,5,0]/netV)**0.5, label='$v_\phi/v$')
	axes[i-1].scatter(dataSelected[:,1], dataSelected[:,8]/dataSelected[:,6], c=a1.get_color())

	axes[i-1].text(0.02,0.05,'$\mathrm{Co}$ = '+str(round(rots[i-1], 3)), transform=axes[i-1].transAxes)
	axes[i-1].grid(False)

	if i < len(rots) - 1:
		axes[i-1].xaxis.set_visible(False)
	axes[i-1].set_ylim([-0.07,1.07])
	axes[i-1].set_xlim([-100,100])
	axes[i-1].locator_params(axis='y', nbins=4)
	axes[i-1].spines['top'].set_visible(False)

axes[-1].set_xlabel('Latitude (degrees)')
axes[-2].set_xlabel('Latitude (degrees)')

h, l = axes[-1].get_legend_handles_labels()
legend = fig.legend(h, l, ncol=3, loc='lower center')

fig.subplots_adjust(wspace=0.12, hspace=0, bottom=0.09)
plt.savefig('Plots/DiagonalComparison.pdf',bbox_inches='tight')

# The third plot compares <v_i v_j> / <v^2> for i != j between our model and both simulations.
# This makes use of Kapyla's Table 3 and Chan's Table 3.

data = np.array(list(KapylaTable3) + list(ChanTable3))

rots = list(set(data[:,0]))
rots.sort()

fig, axes = plt.subplots(len(rots), 1, figsize=(5,12))

latran = np.linspace(-90, 90, num=15, endpoint=True)

for i in range(1,len(rots)+1):
	dataSelected = data[np.abs(data[:,0] - rots[i-1]) < 1e-3]

	res = np.array([f(latitude, rots[i-1]) for latitude in latran])
	netV = res[:,3,3,0] + res[:,4,4,0] + res[:,5,5,0]

	a1, = axes[i-1].plot(latran, res[:,4,3,0]/netV, label='$v_r v_\\theta / v^2 $')
	axes[i-1].scatter(dataSelected[:,1], dataSelected[:,4], c=a1.get_color())

	a1, = axes[i-1].plot(latran, res[:,5,3,0]/netV, label='$v_r v_\phi / v^2 $')
	axes[i-1].scatter(dataSelected[:,1], dataSelected[:,3], c=a1.get_color())

	a1, = axes[i-1].plot(latran, res[:,4,5,0]/netV, label='$v_\\theta v_\phi / v^2 $')
	axes[i-1].scatter(dataSelected[:,1], dataSelected[:,2], c=a1.get_color())

	axes[i-1].text(0.02,0.05,'$\mathrm{Co}$ = '+str(round(rots[i-1], 3)), transform=axes[i-1].transAxes)
	axes[i-1].grid(False)

	axes[i-1].set_xlim([-100,100])
	axes[i-1].set_ylim([-0.23,0.23])

	axes[i-1].locator_params(axis='y', nbins=5)
	axes[i-1].spines['top'].set_visible(False)

	if i == 1:
		axes[i-1].legend(loc='upper left')
	if i < len(rots):
		axes[i-1].spines['bottom'].set_visible(False)
		axes[i-1].xaxis.set_visible(False)

axes[len(rots) - 1].set_xlabel('Latitude (degrees)')
fig.subplots_adjust(wspace=0, hspace=0.0  )
plt.savefig('Plots/OffDiagonalComparison.pdf',bbox_inches='tight')


