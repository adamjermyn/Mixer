import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt
plt.style.use('ggplot')

def paperOne(omega):
	# x is |R grad ln Omega|
	return (6/np.pi)*(omega*(-0.589))*np.sin(np.pi/4)

fi = h5py.File('Data/scale_results.dat','r')
omega = np.array(fi['omega'])
r = np.array(fi['results'])
fi.close()

plt.figure(figsize=(4,3.5))

x, y = omega, r[...,3,5,0]
y = np.abs(y)
plt.loglog(x,y,label='$v_r v_\phi$')
x, y = omega, r[...,3,2,0]
y = np.abs(y)
plt.loglog(x,y,label='$v_r r_\phi$')
plt.loglog(x[x<1], np.abs(paperOne(x[x<1])), label='$v_r v_\phi$ (*)')
plt.xlabel('$ \Omega/|N|$')
plt.ylabel('(Mixing Units)')
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig('Plots/Plot3.pdf')
