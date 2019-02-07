import  os
import  math 
import  numpy              as      np

from    growth_rate        import  growth_factor
from    bz                 import  get_dropoutbz
from    pmm                import  Pmm


def unit_bz(z):
  return  np.ones_like(z)

def linz_bz(z):
  return  (1. + z) 

def invG_bz(z):
  return  1. / growth_factor(z)

def dropout_bz(z):
  bz = get_dropoutbz()
  
  return  bz(z)

def get_PkInterps(cambx):
  Pk_interps              = {}

  Pk_interps['camb_lin']  = cambx.get_lpower_interpolatekz()     ## Interpolates in both z and k.
  Pk_interps['camb_nlin'] = cambx.get_nlpower_interpolatekz()

  return  Pk_interps

def Pmm(Pk_interps, k, z, type = 'nlinear'):
  ##  Required to take an array of k values, and float argument for z.                                                                           
  if   type == 'linear':  
    result   =  Pk_interps['camb_lin'](z, k)

  elif type == 'nlinear':
    result   =  Pk_interps['camb_nlin'](z, k)

  else:
    raise  ValueError('Type %s not available for Pmm.' % type)

  return  result

def Pmh(Pk_interps, bz, k, z):
  return  bz(z) * Pk_interps['camb_nlin'](z, k)

def Phh(Pk_interps, bz, k, z, bz2=None):
  if bz2 is None:
    return  Pk_interps['camb_nlin'](z, k) * bz(z) ** 2.

  else:
    return  Pk_interps['camb_nlin'](z, k) * bz(z) * bz2(z)


if __name__ == "__main__":
  import  pylab      as      pl
  from    prep_camb  import  CAMB


  print("\n\nWelcome.")

  cambx        = CAMB()                                  

  Pk_interps   = get_PkInterps(cambx)

  ks           = np.arange(1e-3, 1.0, 1e-3)

  for z in [2.0, 3.0, 4.0, 5.0]:
    Ps         = Phh(Pk_interps, dropout_bz, ks, z, bz2 = dropout_bz)
    
    pl.loglog(ks, Ps, label = str(z))

  pl.legend()

  pl.show()

  print("\n\nDone.")
