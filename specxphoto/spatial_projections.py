import numpy              as     np
import pylab              as     pl
import astropy.constants  as     const
import matplotlib.pyplot  as     plt

from   scipy.special      import j0, spherical_jn
from   params             import params
from   astropy.cosmology  import FlatLambdaCDM
from   Cgg                import survey_pz, sliced_pz
from   main               import Pab


cosmo            = FlatLambdaCDM(H0 = 100.*params['h_100'], Om0 = params['om_m'], Ob0 = params['om_b'])

def projected_corrfn(Rs, nowiggle=False, alpha=1.0):
  '''                                                                                                                                                                                                                                                     Projected correlation function.                                                                                                                                                                                                                         '''
  dk         = 0.001
  ks         = np.arange(dk,  10.0, dk)

  gks, gRs   = np.meshgrid(ks, Rs)            ## 2D Grids of ks and Rs.                                                                                                                                                                                                  
  Cksp       = projected_crosspower(ks, "SDSS9", "CMASS", nowiggle=nowiggle, alpha=alpha)                                                                                                                                                                                                                              
  ## Fast Bessel function of the first kind of order zero.                                                                                                                                                                                          
  ## integrand  = spherical_jn(0, gks * gRs)
  integrand  = j0(gks * gRs)

  integrand *= ks * Cksp * dk / (2. * np.pi)

  return np.sum(integrand, axis=1)

def projected_crosspower(ks, surveya, surveyb, nowiggle=False, alpha=1.0):
  dz         = 0.01
  zs         = np.arange(dz, 6.0, dz)
  
  gzs, gks   = np.meshgrid(zs, ks)

  prefactor  = (cosmo.H(zs).value/const.c.to('km/s').value) * sliced_pz(zs, 0.75, 1.0, surveya, normed = True) * sliced_pz(zs, 0.75, 0.15, surveyb, normed = True) * dz
  
  integrand  = Pab(gks, ba=1.7, bb=2.0, alpha = alpha, sigma = 0.0, nowiggle=nowiggle)
  integrand *= prefactor

  return np.sum(integrand, axis=1)

def plot_corrfn():
  dR  = 1.000
  Rs  = np.arange(0.001, 145.0, dR)

  for nowiggle in [False]: ## [False, True]:
    for alpha in np.arange(0.90, 1.10, 0.05):
      wRs = projected_corrfn(Rs, nowiggle=nowiggle, alpha=alpha) 

      pl.plot(Rs, Rs * wRs, '-', label="%.2lf" % alpha)

  pl.xlim(0.0,   145.)
  ## pl.ylim(-0.05, 0.30)
  
  pl.xlabel(r'$R \ [h^{-1} \rm{Mpc}]$')
  pl.ylabel(r'$R w(R) \ [h^{-1} \rm{Mpc}]$')

  pl.legend(loc=1, frameon=False)

  pl.savefig('projected_corrfn.pdf', bbox_inches='tight')

def plot_crosspower():
  dk    = 0.0001
  ks    = np.arange(dk, 10.0, dk)
  
  Cksp  = projected_crosspower(ks, "LSST", "CMASS")
  
  pl.loglog(ks, Cksp)

  data  = np.loadtxt('wiggle.dat') 
  pl.plot(data[:,0], data[:,1])

  pl.savefig('projected_crosspower.pdf', bbox_inches='tight')


if __name__ == "__main__":
  ## plot_crosspower()

  plot_corrfn()
