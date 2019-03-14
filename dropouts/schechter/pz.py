import  numpy as np
import  get_wschechters

from    dVols                 import  dVols
from    cosmo                 import  cosmo
from    params                import  get_params
from    goldrush              import  completeness    as grush_completeness
from    Malkan                import  completeness    as malkan_completeness 
from    nbar                  import  dndz
from    get_schechters        import  get_schechters
from    get_wschechters       import  get_wschechter


params = get_params()

def get_pz(zs, ns, C):
  assert  np.allclose(np.diff(zs), np.roll(np.diff(zs), 1))

  zs, Vs  = dVols(zs, cosmo, params, tvol=False)
  Cs      = C(zs)

  ps      = ns[:-1] * Vs * Cs

  dz      = np.unique(np.diff(zs))[0]                                                                                                       

  norm    = np.sum(ps) * dz                                                                                                           
  ps     /= norm                                                                                                                 
  ps     /= dz

  return  zs, ps


if __name__ == '__main__':
  import  pylab as pl

  from    reddy.specs     import  samplestats as reddy_stats
  from    goldrush.specs  import  samplestats as grush_stats
  from    Malkan.specs    import  samplestats as malkan_stats
  from    reddy.pz        import  get_pz      as reddy_getpz


  print('\n\nWelcome to get_schechters.\n\n')

  mlim           =  26.5
  key            =  'g'

  stats          =  grush_stats(printit = True)

  midz, alpha, M_star, phi_star = get_schechters(stats, key)
  ##  midz, alpha, M_star, phi_star     = get_wschechter(4.0)

  zs             =  np.arange(0.0, 10.0, 0.01)
  zs, dVs, ns    =  dndz(zs, phi_star, M_star, alpha, mlim, type='app', printit = True, completeness=None, app_linelim=False)
  
  C              =  grush_completeness.get_completeness(key)

  zs, ps         =  get_pz(zs, ns, C)
  
  ##  Reddy
  ##  zs, ps     =  reddy_getpz(interp=False)

  pl.plot(zs, ps)
  ##  pl.xlim(4.0, 6.0)
  pl.xlabel(r'$z$')
  pl.ylabel(r'$p(z)$')
  pl.show()
  
  print('\n\nDone.\n\n')
