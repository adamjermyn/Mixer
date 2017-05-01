import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt

fi = h5py.File('Data/baro_results.dat','r')
delta = np.array(fi['delta'])
r = np.array(fi['results'])
fi.close()

x, y = delta, r[...,4,3,0]
plt.plot(x,y)
plt.title('Mid-latitude $r \\theta$ stress, $\Omega=0.1|N|$')
plt.xlabel('$\delta$')
plt.ylabel('$\langle v_r v_\\theta \\rangle$')
plt.tight_layout()
plt.savefig('Plots/9_vrvt_baro.pdf')
