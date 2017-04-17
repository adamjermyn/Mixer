import os
from os.path import dirname, abspath
d = dirname(dirname(abspath(__file__)))
os.chdir(d)

import numpy as np
import h5py
import matplotlib.pyplot as plt

fi = h5py.File('Data/scale_results.dat','r')
omega = np.array(fi['omega'])
r = np.array(fi['results'])
fi.close()

plt.subplot(221)

x, y = omega[omega<1], r[omega<1][...,4,4,0]
plt.plot(x,y,label='model')
fit = np.polyfit(x[x < 0.3],y[x < 0.3],2)
plt.plot(x, fit[0]*x**2 + fit[1]*x + fit[2], label='fit: ' + str(round(fit[0],2)) + '*x^2 + ' + str(round(fit[1],2)) + 'x + ' + str(round(fit[2],2)))
plt.legend(loc='upper right')
plt.title('vrvr')
plt.xlabel('$\Omega$')
plt.ylabel('$v_r v_r$')


plt.subplot(222)
x, y = omega[omega<1], r[omega<1][...,4,5,0]
x, y = np.log10(x), np.log10(np.abs(y))
plt.plot(x,y,label='model')
fit = np.polyfit(x[(x > -2.5) & (x < -1)],y[(x > -2.5) & (x < -1)],1)
plt.plot(x, fit[0]*x + fit[1], label='fit: ' + str(round(fit[0],2)) + 'x + ' + str(round(fit[1],2)))
plt.legend(loc='upper right')
plt.title('vrvt')
plt.xlabel('$\log \Omega$')
plt.ylabel('$\log v_r v_\\theta$')


plt.subplot(223)
x, y = omega[omega<1], r[omega<1][...,4,0,0]
plt.plot(x,y,label='model')
fit = np.polyfit(x[x < 0.3],y[x < 0.3],2)
plt.plot(x, fit[0]*x**2 + fit[1]*x + fit[2], label='fit: ' + str(round(fit[0],2)) + '*x^2 + ' + str(round(fit[1],2)) + 'x + ' + str(round(fit[2],2)))
plt.legend(loc='upper right')
plt.title('vrr')
plt.xlabel('$\Omega$')
plt.ylabel('$v_r r$')


plt.subplot(224)
x, y = omega, r[...,4,1,0]
x, y = np.log10(x), np.log10(np.abs(y))
plt.plot(x,y,label='model')
fit = np.polyfit(x[(x > -2.5) & (x < -1)],y[(x > -2.5) & (x < -1)],1)
plt.plot(x, fit[0]*x + fit[1], label='fit: ' + str(round(fit[0],2)) + 'x + ' + str(round(fit[1],2)))
plt.legend(loc='upper right')
plt.title('vrt')
plt.xlabel('$\log \Omega$')
plt.ylabel('$\log v_r \\theta$')

plt.savefig('Plots/RRRT.pdf')
plt.clf()

