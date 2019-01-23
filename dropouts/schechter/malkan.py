import  os
import  numpy             as      np 

from    nbar              import  comovdensity, projdensity
from    Malkan.specs      import  samplestats


def dropout_arealdensity(z, dz, stats, mlims):
  pnbars    = []
  
  alpha     = stats[survey]['schechter']['alpha']
  Mstar     = stats[survey]['schechter']['M_star']
  phi_star  = stats[survey]['schechter']['phi_star']

  print("\n*** Areal density of %s u-band dropout galaxies at z~3. ***"    % survey)
  print("*** Phi_star:  %.6lf;  M_star:  %.6lf;  alpha:  %.6lf ***\n"      % (phi_star, Mstar, alpha))

  for mlim in mlims:
      ## Exluding completeness correction.
      ## projdensity(stats[band]['z'] - dz/2., stats[band]['z'] + dz/2., stats[band]['schechter']['phi_star'], stats[band]['schechter']['M_star'],\
      ## stats[band]['schechter']['alpha'], mlim=mlim, printit = True, completeness=None)  
      
      ## Included HSC completeness.
      pnbar = projdensity(z - dz/2., z + dz/2., phi_star, Mstar, alpha, mlim=mlim,\
                          printit = True, completeness=None)
  
      pnbars.append(pnbar)
      
  np.savetxt('dat/schechter_estimate_%s_dropouts.txt' % survey, np.c_[mlims, np.array(pnbars)], fmt='%.6le')


if __name__ == "__main__":
  print("\n\nWelcome to a Schechter fn. calculator for the areal density of u-band dropouts (Malkan ++). \n\n")

  z         = 3.0
  dz        = 0.4
  
  stats     = samplestats()
  mlims     = np.arange(22., 26.5, 0.1)

  alpha     = stats['Malkan']['schechter']['alpha']
  Mstar     = stats['Malkan']['schechter']['M_star']
  phi_star  = stats['Malkan']['schechter']['phi_star']
  '''
  for survey in ['Malkan']:
    ## dropout_arealdensity(z, dz, stats, mlims)

    for mlim in np.arange(23.0, 26.0, 0.1):
      comovdensity(z, phi_star, Mstar, alpha, type='app', mlim=mlim, printit=True)
  '''
  print("\n\nDone.\n\n")
