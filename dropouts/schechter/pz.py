import  numpy as np

from    dVols                 import  dVols
from    cosmo                 import  cosmo
from    params                import  get_params
from    goldrush              import  completeness 
from    nbar                  import  dndz


params = get_params()

def get_pz(zs, ns, C):
  ##  Note:  returns normalised completeness curve.                                                                                                      
  ##         requires LF and volume factors to be added.                                                                                                 
  params         = get_params()

  zs, Vs         = dVols(zs, cosmo, params, tvol=False)

  Cs             = C(zs)

  I              = ns * Vs * Cs

  ##  dz         = zee[1] - zee[0]                                                                                                                       
  ##  norm       = np.sum(completeness * dz)                                                                                                               

  ##  pzee       = completeness / norm                                                                                                                    

  ##  return  zee, pzee    


if __name__ == '__main__':
  zs             = np.arange(0.0, 10.0, 0.01)
  zs, dVs, ns    = dndz(zs, phi_star, M_star, alpha, mlim, type='app', printit = True, completeness=None, app_linelim=False)

  C              = completeness.get_completeness(drop)

  ##  get_pz()

  print('\n\nDone.\n\n')
