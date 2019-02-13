import numpy as np
import pylab as pl


dat   = np.loadtxt('rsd.dat')
dat   = dat[1::2]

k1    = dat[::2]
k2    = dat[1::2]

zmean = (k1[:,1] + k1[:,2]) / 2.

print(k1)
print(k2)

pl.plot(zmean[zmean < 2.0], k1[:,4][zmean < 2.0], c='r', label=r'$k_{\rm{max}}=0.1$ [$h$/Mpc]', marker='^')
pl.plot(zmean[zmean < 2.0], k2[:,4][zmean < 2.0], c='r', label=r'$k_{\rm{max}}=0.2$ [$h$/Mpc]', marker='s')

pl.plot(zmean[zmean > 2.0], k1[:,4][zmean > 2.0], c='orange', label=r'$k_{\rm{max}}=0.1$ [$h$/Mpc]', marker='^')
pl.plot(zmean[zmean > 2.0], k2[:,4][zmean > 2.0], c='orange', label=r'$k_{\rm{max}}=0.2$ [$h$/Mpc]', marker='s')

pl.xlabel(r'$z$')
pl.ylabel(r'$\sigma_f$')

pl.legend(ncol=2)
pl.savefig('sigf.pdf')
