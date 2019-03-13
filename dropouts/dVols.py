import numpy  as     np

from   cosmo  import cosmo
from   params import get_params


params = get_params()

def dVols(zs, cosmo, params, tvol=False):
  Vs    = cosmo.comoving_volume(zs).value                                           ##  Get volume to each redshift slice.                                
  Vs   *= params['h_100']**3.                                                       ##  [h^-1 Mpc]^3                                                      

  dVs   = Vs - np.roll(Vs, 1)

  ##  Lower z limit on each slice.                                                                                                                      
  zs    = zs[:-1]
  dVs   = dVs[1:]

  if tvol:
    ##  Total volume of the slice.                                                                                                                       
    return  zs, dVs, Vs[-1] - Vs[0]

  else:
    return  zs, dVs


if __name__ == '__main__':
    import pylab as pl


    zs     = np.arange(0.0, 10.0, 0.1)
    zs, Vs = dVols(zs, cosmo, params, tvol=False)
    
    pl.plot(zs, Vs)
    pl.show()

    print('\n\nDone.\n\n')
