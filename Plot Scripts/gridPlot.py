import numpy as np
import h5py
import matplotlib.pyplot as plt

fi = h5py.File('../Data/results.dat','r')
omega = np.array(fi[u'omega'])
omega = np.log10(omega)
theta = np.array(fi[u'theta'])
wOverOmega = np.array(fi[u'wOverOmega'])
wOverOmega = np.log10(wOverOmega)
results = np.array(fi[u'results'])
fi.close()

vrvr = results[...,4,4,0]
vrvtheta = results[...,5,4,0]
vrvp = results[...,6,4,0]
vrr = results[...,4,0,0]
vrt = results[...,4,1,0]

# 1
plt.plot(10**omega, vrvr[:,len(theta) / 4,0])
plt.title('Mid-Latitude $r r$ Stress')
plt.xlabel('$ \Omega/|N|$')
plt.ylabel('$\langle v_r v_r \\rangle$')
plt.savefig('Plots/vrvr_omega_1.pdf')
plt.clf()

# 2
plt.plot(omega, np.log10(np.abs(vrvtheta[:,len(theta) / 4,0])))
plt.title('Mid-latitude $r \\theta$ stress')
plt.xlabel('$\log \Omega/|N|$')
plt.ylabel('$\log |\langle v_r v_\\theta\\rangle|$')
plt.savefig('Plots/vrvt_omega_2.pdf')
plt.clf()

# 3
plt.plot(omega, np.log10(np.abs(vrvp[:,len(theta) / 4,0])))
plt.title('Mid-latitude $r \phi$ stress')
plt.xlabel('$\log \Omega/|N|$')
plt.ylabel('$\log|\langle v_r v_\phi\\rangle|$')
plt.savefig('Plots/vrvp_omega_3.pdf')
plt.clf()

# 4
plt.plot(omega, np.log10(np.abs(vrvp[len(omega) / 2,len(theta) / 4,:])))
plt.title('Mid-latitude $r \phi$ stress')
plt.xlabel('$\log |\\nabla \log \Omega |$')
plt.ylabel('$\log|\langle v_r v_\phi \\rangle|$')
plt.savefig('Plots/vrvp_grad_omega_4.pdf')
plt.clf()

# 5
plt.plot(omega, np.log10(np.abs(vrvtheta[len(omega) / 2,len(theta) / 4,:])))
plt.title('Mid-latitude $r \\theta$ stress')
plt.xlabel('$\log |\\nabla \log \Omega |$')
plt.ylabel('$\log|\langle v_r v_{\\theta} \\rangle|$')
plt.savefig('Plots/vrvt_grad_omega_5.pdf')
plt.clf()

# 6
plt.plot(10**omega, vrr[:,len(theta) / 4,0])
plt.title('Mid-latitude $r r$ viscosity')
plt.xlabel('$\Omega/|N|$')
plt.ylabel('$\langle v_r \delta r \\rangle$')
plt.savefig('Plots/vrr_omega_6.pdf')
plt.clf()

# 7
plt.plot(omega, np.log10(np.abs(vrt[:,len(theta) / 4,0])))
plt.title('Mid-latitude $r \\theta$ viscosity')
plt.xlabel('$\log \Omega/|N|$')
plt.ylabel('$\log|\langle v_r \delta (r \\theta) \\rangle|$')
plt.savefig('Plots/vrt_omega_7.pdf')
plt.clf()




