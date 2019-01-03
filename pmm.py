import  os
import  math 
import  numpy              as      np

from    growth_rate        import  growth_factor


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
    result   = Pk_interps['camb_nlin'](z, k)

  else:
    raise  ValueError('Type %s not available for Pmm.' % type)

  return  result


if __name__ == "__main__":
  import  pylab      as      pl
  from    prep_camb  import  CAMB


  print "\n\nWelcome."

  cambx        = CAMB()                                  

  Pk_interps   = get_PkInterps(cambx)

  ks           = np.arange(1e-3, 10.0, 1e-3)
  Ps           = Pmm(Pk_interps, ks, 0.6, type = 'linear')

  pl.loglog(ks, Ps)
  pl.show()

  np.savetxt('dat/lin_pmm.dat', np.c_[ks, Ps], fmt='%.6le')

  print("\n\nDone.")
