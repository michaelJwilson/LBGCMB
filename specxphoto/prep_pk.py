import  numpy as np

from    growth_rate        import  growth_factor
from    prep_camb          import  CAMB
from    pmh                import  Pmm, get_PkInterps, linz_bz


def prep_pk(z):
  try:
    ks, Ps = np.loadtxt('dat/pk.dat', unpack=True)

  except:
    ##  Prepare pycamb module; linear, non-linear matter P(k) and Cls.                                                                                     
    cambx       =  CAMB()
    Pk_interps  =  get_PkInterps(cambx)

    ##  Write Linear power spectrum.                                                                                                                     
    ks = np.arange(0.01, 1., 0.01)
    Ps = Pk_interps['camb_lin'](0.0, ks)

    np.savetxt('dat/pk.dat', np.c_[ks, Ps])

  aa = 1. / (1. + z)
  dd = growth_factor(aa)

  return  ks, dd * dd * Ps
