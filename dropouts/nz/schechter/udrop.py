import  numpy  as      np
import  pylab  as      pl

from    nbar   import  comovdensity, projdensity


if __name__ == "__main__":
  from  reddy  import  samplestats

  """                                                                                                                                                      
  Reddy estimates will be biased as sample was targeted on QSOs,                                                                                           
  i.e. volumes of known over-density.                                                                                                                      
  """

  print("\n\nWelcome to the Schechter fn. estimate of Reddy u-drops.")


  stats       =  samplestats(printit = True)

  z           =  stats['LBG']['z']
  dz          =  0.5

  alpha       =  stats['LBG']['alpha']
  M_star      =  stats['LBG']['M_star']
  phi_star    =  stats['LBG']['phi_star']

  ## ---------------------------------------------------------------------- ##
  print "\n\n*** Reddy & Steidel z~3 u-dropout survey ***"
  print "*** Comoving number density of u-band dropout galaxies at z~3. ***"
  print "\nR_lim            M_UV                    L/L*            log10(nbar[(h_100/Mpc)^3])"

  for Rlim in np.arange(23., 26.5, 0.5):
    comovdensity(z, phi_star, M_star, alpha, type='app', mlim=Rlim, printit = True)

  ## ---------------------------------------------------------------------- ##
  print "\n\n*** Areal density of u-band dropout galaxies at z~3. ***"
  print "\nR_lim                            Galaxies / deg.^2               Volume / Gpc^3"

  for Rlim in np.arange(23., 26.5, 0.5):
    projdensity(z-dz/2., z+dz/2., phi_star, M_star, alpha, mlim=Rlim, printit = True)
    
  print("\n\nDone.\n\n")
