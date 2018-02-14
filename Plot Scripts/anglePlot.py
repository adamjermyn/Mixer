import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

from KRflux import flux

import numpy as np
import h5py
import matplotlib.pyplot as plt
plt.style.use('ggplot')

fi = h5py.File('Data/scale_results.dat','r')
omega = np.array(fi['omega'])
r = np.array(fi['results'])
fi.close()

print(r[4,:,:,0])

plt.figure(figsize=(6,4.3))

x, y, z = omega, r[...,0,3,0], r[...,3,4,0]
plt.plot(x,np.arctan2(z,y),label='$\\theta$')
plt.xscale('symlog',linthreshx=1.0)
plt.ylim([-0.21,0.02])
plt.xlabel('$ \Omega/|N|$')
plt.ylabel('Baroclinic Angle (Rad)')
plt.tight_layout()
plt.savefig('Plots/AnglePlot.pdf')
