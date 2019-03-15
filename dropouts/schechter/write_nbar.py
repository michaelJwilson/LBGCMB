import  os
import  pylab                                         as pl
import  numpy                                         as np
import  get_wschechters

from    dVols                 import  dVols
from    cosmo                 import  cosmo
from    params                import  get_params
from    goldrush              import  completeness    as grush_completeness
from    Malkan                import  completeness    as malkan_completeness
from    nbar                  import  dndz
from    get_schechters        import  get_schechters
from    get_wschechters       import  get_wschechter
from    scipy.interpolate     import  interp1d
from    reddy.specs           import  samplestats     as reddy_stats
from    goldrush.specs        import  samplestats     as grush_stats
from    Malkan.specs          import  samplestats     as malkan_stats
from    reddy.pz              import  get_pz          as reddy_getpz
from    nbar                  import  comovdensity


if __name__ == '__main__':
  print('\n\nWelcome: write number densities for each sample with mag. limit.')

  band                            =   'u' 
  specfrac                        =  {'BX': 1.0, 'u': 1.0, 'g': 0.1, 'r': 0.1}

  if band == 'BX':
    stats                         =   reddy_stats(printit = False)
    midz, alpha, M_star, phi_star = get_schechters(stats, band)

  elif band == 'u':
    stats                         = malkan_stats(printit  = False)    
    midz, alpha, M_star, phi_star = get_schechters(stats, 'Malkan')

  else:
    stats                         =  grush_stats(printit = False)
    midz, alpha, M_star, phi_star = get_schechters(stats,  band)  ##  midz, alpha, M_star, phi_star = get_wschechter(4.0)                                    

  ## 
  mlims                           = np.arange(22.5, 30.0, 0.05)
  result                          = []

  for mlim in mlims:
    result.append(comovdensity(midz, phi_star, M_star, alpha, type='app', mlim=mlim, printit=False))
  
  result                          = np.array(result)

  np.savetxt(os.environ['LBGCMB'] + '/dropouts/schechter/dat/comov_density/%sDrop.dat' % band, np.c_[mlims, result, np.log10(specfrac[band]) + result],\
                                                                                               fmt='%.6lf \t %.6le \t %.6le', header='mlim  photo  spec')

  print('\n\nDone.\n\n')
