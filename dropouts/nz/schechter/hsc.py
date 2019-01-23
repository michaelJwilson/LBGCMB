import  os
import  numpy           as      np 

from    nbar            import  comovdensity, projdensity
from    goldrush.specs  import  samplestats
from    completeness    import  interp_completeness
from    goldrush.ilim   import  get_nbar_nocontam


def dropout_arealdensity(z, dz, band, stats, mlims):
  alpha     =  stats[band]['schechter']['alpha']
  Mstar     =  stats[band]['schechter']['M_star']
  phi_star  =  stats[band]['schechter']['phi_star']

  pnbars    = []

  print("\n*** Areal density of HSC %s-band dropout galaxies at z~%d. ***" % (band, np.round(stats[band]['z'])))
  print("*** Phi_star:  %.6lf;  M_star:  %.6lf;  alpha:  %.6lf ***\n"      % (phi_star, Mstar, alpha))
  
  for mlim in mlims:
      ##  Exluding completeness correction.
      ##  projdensity(stats[band]['z'] - dz/2., stats[band]['z'] + dz/2., stats[band]['schechter']['phi_star'], stats[band]['schechter']['M_star'],\
      ##  stats[band]['schechter']['alpha'], mlim=mlim, printit = True, completeness=None)  
      
      ##  Included HSC completeness.
      pnbar = projdensity(stats[band]['z'] - dz / 2., stats[band]['z'] + dz / 2., phi_star, Mstar, alpha, mlim=mlim,\
                          printit=True, completeness=interp_completeness)  

      pnbars.append(pnbar)

  np.savetxt('dat/schechter_estimate_%s_dropouts.txt' % band, np.c_[mlims, np.array(pnbars)], fmt='%.6le')

  stats = get_nbar_nocontam(band, depth='W', printit=False)   
  
  print("\nGoldRush observed:  %.3le g/deg2; Contamination corrected:  %.3le g/deg2" % (stats[band]['nbar'], stats[band]['nbar_nointerlopers']))


if __name__ == "__main__":
  print("\n\nWelcome to a Schechter fn. calculator for the areal density of HSC Goldrush dropouts.\n\n")
  
  band      = 'g'
  stats     = samplestats()
  
  z         = stats[band]['z']

  ##  Redshift distribution limits handles by completeness(z), rather than dz.
  dz        = 1.5
  
  alpha     = stats[band]['schechter']['alpha']
  Mstar     = stats[band]['schechter']['M_star']
  phi_star  = stats[band]['schechter']['phi_star']

  mlims     = np.arange(22., 26.5, 0.1)

  for band in ['g']:
    dropout_arealdensity(z, dz, band, stats, mlims)

    ##  for mlim in np.arange(23.0, 26.0, 0.1):
    ##    comovdensity(z, phi_star, Mstar, alpha, type='app', mlim=mlim, printit=True)

  print("\n\nDone.\n\n")
