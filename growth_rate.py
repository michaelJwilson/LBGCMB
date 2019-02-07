import numpy   as      np

from   cosmo   import  cosmo
from   params  import  get_params


params = get_params()

def growth_rate(a): 
  gamma = 0.545
  z     = 1./a - 1.

  # rho_c   = 3 * H^2 / (8 * pi G)                                                                                                                         
  # om_m(a) = (rho_m/rho_c) = 8 * pi * G *rho_m(0) * a^-3 / (3 H^2)                                                                                        
  #                         = om_m(a = 1) * a^-3 * (H_0/H)^2                                                                                                

  om_ma = params['om_m'] * (cosmo.H(0) / cosmo.H(z))**2. * a**-3.

  return  om_ma.value**gamma

def _growth_factor(aa):  
  avals     =  np.arange(0.001, 1., 0.001)
  integrand = (growth_rate(avals) - 1.) / avals

  norm      =    np.exp(np.trapz(integrand,              avals))
  result    = aa*np.exp(np.trapz(integrand[avals <= aa], avals[avals <= aa]))
  
  return result / norm

growth_factor = np.vectorize(_growth_factor)

def plot_dplus():
  import pylab as pl

  pl.clf()
      
  zs = np.arange(0.001, 1., 0.01)
  aa = 1./(1. + zs)

  pl.plot(zs,  growth_factor(aa), 'k-')
  pl.plot(zs,                aa,  'r--')

  pl.xlabel('z')
  
  pl.savefig('plots/dplus.pdf')


if __name__ == "__main__":
  plot_dplus()

  ## print 0.84 / _growth_factor(1. / (1. + 1.13))
