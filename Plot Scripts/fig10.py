import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt

fi = h5py.File('Data/accretion_results.dat','r')
B = np.array(fi['B'])
r = np.array(fi['results'])
fi.close()

x, y = np.log10(B), r[...,3,5,0]
plt.plot(x,y)
plt.title('Mid-latitude $r \phi$ stress $|N|=0$')
plt.xlabel('$\log v_A/l\Omega$')
plt.ylabel('$\langle v_r v_\phi \\rangle$')
plt.tight_layout()
plt.savefig('Plots/10_vrvp_accrete.pdf')
