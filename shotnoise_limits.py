import numpy as np
import pylab as pl


def shotnoise_limit(Llls, Cgg, fsky = 0.5, lmax = 1500., unit='arcmin', printit = True, type = 'equal', redshift=None):
  if type == 'equal':
    Ngg   = Cgg

  else:
    raise ValueError("Type is not available for shotnoise limit.")

  nbar    = 1./Ngg                          ## galaxies per steradian.                                                                               

  deg_str = np.pi/180.

  nbar   *= deg_str**2.                     ## galaxies per sq. deg.                                                                                                                                                                                                                                                        
  Index   = (np.abs(Llls - lmax)).argmin()

  nbar    = nbar[Index]                     ## shot noise limit @ lmax.                                                                                                                                                                                                                                              
  if printit == True and type == 'equal':
    print "\n\nFor Ngg to equal Cgg @ Lmax of %d (z of %.3lf) requires %.3lf (%.3lf) galaxies per sq. deg. (arcmin)\n"   % (lmax, redshift, nbar, nbar/60./60.)

  return nbar

def signal2noise(Llls, Ckg, Ckk, Nkk, Cgg, Ngg, fsky = 0.5, lmax = 1500.):
  ''' Return signal to noise as a function of Ngg. '''

  r_ell  = Ckg / np.sqrt((Ckk + Nkk) * (Cgg + Ngg))

  s2n    = 1./np.sqrt(1. + r_ell ** -2.)   ## per ell                                                                                                    
  s2n   *=  np.sqrt(fsky*(2.*Llls + 1.))   ## Accounting for all ms and area.                                                                  
  
  Index  = (np.abs(Llls - lmax)).argmin()  ## L closest to lmax.                                                                                  
  
  return  s2n[Index]


if __name__ == "__main__":
  print "\n\nWelcome."
