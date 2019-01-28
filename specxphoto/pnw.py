import numpy              as      np
import pylab              as      pl

from   qso_bz             import  qso_bz
from   pmh                import  linz_bz
from   scipy.interpolate  import  spline
from   prep_pk            import  prep_pk


def pnw(z, b1, b2):
  ##  No wiggle spectrum from k^-3 to k^2 spline. 
  ks, ps   = prep_pk(z)

  fitz     = np.polyfit(np.log10(ks), np.log10(ps / ks ** 3.), 5)
  ps       = ks ** 3. * 10. ** np.poly1d(fitz)(np.log10(ks))

  return  ks, b1(z) * b2(z) * ps


if __name__ == '__main__':
  from main import Psp


  print("\n\nWelcome to Nishizwa.\n\n")

  zspec      = 1.5
  ks, ps     = prep_pk(zspec)

  ks, pnws   = pnw(zspec, ks)

  pl.semilogx(ks, ps / pnws)

  pl.legend()
  pl.show()

  print('\n\nDone.\n\n')
