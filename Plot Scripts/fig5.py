import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt

fi = h5py.File('Data/scale_grad_results.dat','r')
w = np.array(fi['w'])
omega = np.array(fi['omega'])
r = np.array(fi['results'])
fi.close()

x, y = w/omega, r[...,5,3,0]
plt.plot(x,y)
plt.title('Mid-latitude $r \phi$ stress, $\Omega=0.1|N|$')
plt.xlabel('$R\\nabla \ln \Omega $')
plt.ylabel('$\langle v_r v_\phi \\rangle$')
plt.tight_layout()
plt.savefig('Plots/5_vrvp_grad_omega.pdf')
