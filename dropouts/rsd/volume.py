def vol_integrand(z, fsky=0.5, fkp_weighted=False, nbar=1.e-3, P0=5.e3):
    ##  Purely volume integral as a sanity check.                                                                                                                              
    ##  dV / dz [(h^{-1} Mpc)^3];  Differential comoving volume per redshift per steradian.                                                                                  
    ##                                                                                                                                                                             
    ##  Note:  cosmo.differential_comoving_volume(z) = (const.c.to('km/s') / cosmo.H(z)) * cosmo.comoving_distance(z) ** 2.                                                      
    ##                                               = dV/dz [d\Omega] = chi^2 dChi/dz [d\Omega].                                                                                
    ##                                                                                                                                                                            
    dVdz = fsky * 4. * np.pi * cosmo.differential_comoving_volume(z).value * cparams['h_100'] ** 3.

    if fkp_weighted:
      ##  FKP volume weighting.                                                                                                                                                
      nP  = nbar * P0
      fkp = nP / (1. + nP)

      return  fkp * fkp * dVdz

    else:
      return  dVdz

def _vvol_integrand(x, args):
    ##  Vegas wrapper of vol_integrand; input args as a list.                                                                                                                     
    z                              = x[0]
    (fsky, fkp_weighted, nbar, P0) = args

    return  vol_integrand(z, fsky, fkp_weighted, nbar, P0)
