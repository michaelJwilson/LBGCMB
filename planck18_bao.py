import  matplotlib         as      mpl
import  matplotlib.pyplot  as      plt

import  numpy              as      np
import  pylab              as      pl 

from    pylab              import  rcParams
from    utils              import  latexify
from    scipy.interpolate  import  interp1d


def get_sig8z(interp=False):
  data   = np.loadtxt('./dat/planck18_bao.dat')

  zs     = data[:,0]
  sig8   = data[:,1]
  esig8  = data[:,2]

  if interp:
    return  zs, interp1d(zs, sig8, kind='linear', bounds_error=True, assume_sorted=False), interp1d(zs, esig8, kind='linear', bounds_error=True, assume_sorted=False) 

  else:
    return  zs, sig8, esig8

def get_fsig8(interp=False):
  data    = np.loadtxt('./dat/planck18_bao.dat')

  zs      = data[:,0]
  fsig8   = data[:,3]
  efsig8  = data[:,4]

  if interp:
    return  zs, interp1d(zs, fsig8, kind='linear', bounds_error=True, assume_sorted=False), interp1d(zs, efsig8, kind='linear', bounds_error=True, assume_sorted=False)

  return  zs, fsig8, efsig8  


if __name__ == '__main__':
    print('\n\nWelcome to Planck18 + BAO.\n\n')

    latexify(fig_width=None, fig_height=None, columns=1, equal=True, usetex=True, fontsize=10)

    zs,  sig8,  esig8 = get_sig8z()
    zs, fsig8, efsig8 = get_fsig8()

    pl.errorbar(zs,          sig8,   esig8, fmt='^', label=r'$\sigma_8(z)$',        capsize=5, markersize=3, linestyle='', alpha=0.6)
    pl.errorbar(zs + 0.085, fsig8,  efsig8, fmt='s', label=r'$f(z) \ \sigma_8(z)$', capsize=5, markersize=3, linestyle='', alpha=0.6)

    pl.xlim(1.8, 5.2)
    pl.ylim(.15, .35)

    pl.xlabel(r'$z$')
    pl.legend(loc=1, frameon=False)

    ##  pl.show()
    pl.savefig('plots/planck18_bao.pdf', bbox_inches='tight')

    print('\n\nDone.\n\n')
