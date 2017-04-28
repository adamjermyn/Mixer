import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt

fi = h5py.File('Data/polar_results.dat','r')
B = np.array(fi['B'])
r = np.array(fi['results'])
fi.close()

print(r[...,0,4,4,0])
print(r[...,1,4,4,0])
print(r[...,2,4,4,0])

x, y = np.log10(B), np.log10(np.abs(r[...,0,4,4,0]))
plt.plot(x,y,label='$r$')
x, y = np.log10(B), np.log10(np.abs(r[...,1,4,4,0]))
plt.plot(x,y,label='$\\theta$')
x, y = np.log10(B), np.log10(np.abs(r[...,2,4,4,0]))
plt.plot(x,y,label='$\phi$')
plt.legend()
plt.title('Mid-latitude $r r$ stress verus magnetic polarisation')
plt.xlabel('$\log v_A/l$')
plt.ylabel('$\log |\langle v_r v_r \\rangle |$')
plt.tight_layout()
plt.savefig('Plots/8_vrvr_B.pdf')
