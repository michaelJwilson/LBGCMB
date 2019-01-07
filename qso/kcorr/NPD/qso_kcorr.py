import os 
import numpy as np
import pylab as pl

from   matplotlib.ticker import AutoMinorLocator


def get_qsokcorr(plotit=False):
  ##  Fig. 9 of https://arxiv.org/pdf/1509.05607.pdf
  zs  = 0.1 * np.array([6., 7., 8.0, 9.0, 10.0, 11., 12., 13., 14., 15., 16., 17., 18., 19.,   20.,   21.,  22.,  23., 24.,  25., 26., 27.,   28., 29.,   30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40.0])
  Ks  = 0.1 * np.array([5., 4., 3.5, 3.2,  3.5, 3.2, 3.2, 3.,  2.5,  2.,  1., 0.5, 0.5, 0.0, -0.25,  0.25, -0.5, -1.5, -2., -1.5, -1., -0.8, -0.7, -0.6, -0.5,  0.,  0., 2.,   5.,  6.,  7.,  8.,  9., 11., 12.5])
  Ksp =       np.polyfit(zs, Ks, 12)
  Ksp =       np.poly1d(Ksp)

  if plotit:
    pl.plot(zs, Ks,      'k^', markersize=3, label='Palanque-Delabrouille 2015')
    pl.plot(zs, Ksp(zs), 'r-')

    pl.xlim(0.6, 4.0)
    pl.ylim(-0.7, 1.7)

    pl.xlabel(r'$z$')
    pl.ylabel(r'$K(z)$')
    pl.legend()

    ax = pl.gca()
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    pl.show()

  return zs, Ks, Ksp


if __name__ == '__main__':
  zs, Ks, Ksp = get_qsokcorr(plotit=True)

  root = os.environ['LBGCMB']

  np.savetxt(root + '/qso/kcorr/NPD/NPD-kcorr.txt', np.c_[zs, Ks, Ksp(zs)], fmt='%.6lf \t %.6lf \t %.6lf')
