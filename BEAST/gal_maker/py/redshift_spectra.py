import  numpy  as      np 
from    utils  import  luminosity_distance


def redshift_spectra(vs, Fv, redshift):
  """                                                                                                                                                        
  vs: restframe / emitted; Similarly, Fv [ergs/s/Hz].                                                                                                        
  eqn. (16.3) on pg. 691 of Mo, Bosch and White.                                                                                                             
  """
  
  Dl   =  luminosity_distance(redshift, unit = 'cm')                                   ## Luminosity distance in cm.                              
  Fv   =  Fv / Dl ** 2.

  Fv  /=  (4. * np.pi)
  Fv  *=  (1. + redshift)

  vs   =  vs / (1. + redshift)                                                         ## Redshift to observed wavelengths; ls are rest frame.        
  
  return  vs, Fv
