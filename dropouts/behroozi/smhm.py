import numpy             as      np
import pylab             as      pl

import matplotlib.pyplot as      plt

from   scipy.misc        import  derivative
from   utils             import  latexify
from   params            import  get_params  as get_cparams


cparams = get_cparams()

def _f(x, alpha, delta, gamma):
  num = np.log10(1. + np.exp(x)) ** gamma
  den = 1. + np.exp(10. ** -x)

  return  -np.log10(10. ** (alpha * x) + 1.) + delta * num / den

def _dfdx(x, alpha, delta, gamma):
  ##  https://tinyurl.com/y7a9hsl5
  interim  = - alpha * 10. ** (alpha * x) / (1. + 10. ** (alpha * x))
  interim += delta * 10. ** (-x) * np.exp(1) ** (10. ** (-x)) * np.log(10.) ** (1. - gamma) * np.log(np.exp(x) + 1.) ** gamma / (1. + np.exp(1) ** (10. ** (-x)))**2.
  interim += gamma * delta * np.exp(x) * np.log(10.) ** (-gamma) * np.log(np.exp(x) + 1.) ** (gamma-1) / (1. + np.exp(1) ** (10.**(-x))) / (np.exp(x) + 1.)

  return interim

def _xc(alpha, delta, gamma, plotit=False):
  ##  Find the value of x=np.log10(M/M1) that provides a stationary point for M*/M.
  ##  This is df/dx = 1.
  dx          = 1.e-5
  xs          = np.arange(dx, 1., dx)

  analytic    = _dfdx(xs, alpha, delta, gamma) 
  ngradient   = np.gradient(_f(xs, alpha, delta, gamma), dx)
  sderivative = derivative(_f, xs, dx=dx, n=1, args=(alpha, delta, gamma), order=3)

  if plotit:
    ## pl.semilogy(xs, ngradient, '^', label='numeric', markersize=2)                                                                                                                          
    pl.semilogx(xs, analytic, label='analytic')
    pl.semilogx(xs, sderivative, '--', label='scipy')

  ##  Find closest x to df/dx = 1.
  index = np.where(np.abs(analytic - 1.).min() == np.abs(analytic - 1.))   
  xc    = xs[index]

  return  xc

def _Mc(alpha, delta, gamma, M1):
    ##  Return the Mass that maximises the stellar to 
    ##  halo mass relation.

    xc = _xc(alpha, delta, gamma, plotit=False)

    return M1 * 10. ** xc  

def get_params():
    ##  Secton 5 of Behroozi 2013.     
    params              =     {}

    params['epsilon0']  = -1.777
    params['epsilona']  = -0.006
    params['epsilonz']  =  0.000
    params['epsilona2'] = -0.119

    params['alpha0']    = -1.412
    params['alphaa']    =  0.731

    params['delta0']    =  3.508
    params['deltaa']    =  2.608
    params['deltaz']    = -0.043

    params['gamma0']    =  0.316
    params['gammaa']    =  1.319
    params['gammaz']    =  0.279

    params['MICL0']     =  12.515
    params['MICLa']     =  -2.503

    params['Sigma0']    =  0.2180
    params['Sigmaa']    =  -0.023

    params['M10']       =  11.514
    params['M1a']       = -1.7930
    params['M1z']       = -0.2510

    ##  Stellar mass completeness.
    params['A']         =  0.273
    params['zc']        =  1.077

    return  params

def _params(z):
  a       = 1. / (1. + z)
  nu      = np.exp(- 4. * a * a)
  
  params  = get_params()

  epsilon = params['epsilon0'] + (params['epsilona'] * (a - 1.) + params['epsilonz'] * z) * nu + params['epsilona2'] * (a-1.)
  epsilon = 10. ** epsilon

  alpha   = params['alpha0'] +  params['alphaa'] * (a-1.) * nu
  delta   = params['delta0'] + (params['deltaa'] * (a-1.) + params['deltaz'] * z) * nu
  gamma   = params['gamma0'] + (params['gammaa'] * (a-1.) + params['gammaz'] * z) * nu

  ## 0.2 +- 0.03 dex. 
  Sigma   = params['Sigma0'] + params['Sigmaa'] * (a - 1.)
  
  M1      = params['M10'] + (params['M1a'] * (a-1.) + params['M1z'] * z) * nu
  M1      = 10. ** M1

  ##  Survey completeness.
  ci      = lambda zee: 1.0 - params['A'] / (1. + np.exp(params['zc'] - zee))
  c       = lambda zee: 1.0 if zee < 1.0 else ci(zee) + 1. - ci(1.0)

  return  epsilon, alpha, delta, gamma, Sigma, M1, c(z) 

def ishikawa(drop, mhalo):
  ##  Ishikawa h is 0.7;
  ##  Top right of pg. 2.
  
  if drop == 'u':
    epsilon, alpha, delta, gamma, Sigma, M1, c = _params(3.0)
    xc       = _xc(alpha, delta, gamma, plotit=False)

    ##  Last value is the halo mass that maximises M*/Mh relation.
    print('%.6le, %.6le, %.6le' % (epsilon, M1, M1 * 10. ** xc))

    ##  Rewrite with Isikawa values; Table 3 of 1612.06869. 
    ishh     = 0.7   ## Hubble  
    ishMc    = 12.10
    isheps   = 2.222e-2

    Mc       = (10. ** ishMc) / ishh
    M1       = Mc * 10. ** (-xc) 

    ##  Provides M* = ishepsilon * ishMc | Mc = M1 * 10 ** xc.  Eavlue at this point with eqn. (3) of 1207.6105 allows
    ##  solving for Behroozi epsilon.
    epsilon  = 10. ** (np.log10(isheps * Mc) + _f(0, alpha, delta, gamma) - _f(xc, alpha, delta, gamma)) / M1

    print('%.6le, %.6le, %.6le, %.6le' % (xc, Mc, epsilon, M1))

  elif drop == 'g':
    epsilon, alpha, delta, gamma, Sigma, M1, c = _params(4.0)
    xc       = _xc(alpha, delta, gamma, plotit=False)

    ##  Last value is the halo mass that maximises M*/Mh relation.                                                                                                                                                             
    print('%.6le, %.6le, %.6le' % (epsilon, M1, M1 * 10. ** xc))
    
    ##  Rewrite with Isikawa values; Table 3 of 1612.06869.                                                                                                                                                                    
    ishh     = 0.7   ## Hubble                                                                                                                                                                                                 
    ishMc    = 11.99 
    isheps   = 2.248e-2

    Mc       = (10. ** ishMc) / ishh                                                                                                                                                                          
    M1       = Mc * 10. ** (-xc)

    epsilon  = 10. ** (np.log10(isheps * Mc) + _f(0, alpha, delta, gamma) - _f(xc, alpha, delta, gamma)) / M1
      
    print('%.6le, %.6le, %.6le, %.6le' % (xc, Mc, epsilon, M1))

  elif drop == 'r':
    epsilon, alpha, delta, gamma, Sigma, M1, c = _params(5.0)
    xc       = _xc(alpha, delta, gamma, plotit=False)

    ##  Last value is the halo mass that maximises M*/Mh relation.                                                                                                                                                             
    print('%.6le, %.6le, %.6le' % (epsilon, M1, M1 * 10. ** xc))

    ##  Rewrite with Isikawa values; Table 3 of 1612.06869.                                                                                                                                                                    
    ishh     = 0.7   ## Hubble                                                                                                                                                                                                 
    ishMc    = 11.77
    isheps   = 2.091e-2

    Mc       = (10. ** ishMc) / ishh
    M1       = Mc * 10. ** (-xc)

    epsilon  = 10. ** (np.log10(isheps * Mc) + _f(0, alpha, delta, gamma) - _f(xc, alpha, delta, gamma)) / M1

    print('%.6le, %.6le, %.6le, %.6le' % (xc, Mc, epsilon, M1))

  else:
    raise ValueError('Requested drop type is not available for Isikawa (2017); 1612.06869')
  
  x  = np.log10(mhalo / M1)
  SM = np.log10(epsilon * M1) + _f(x, alpha, delta, gamma) - _f(0.0, alpha, delta, gamma)

  return  10. ** SM, Sigma, c

def _med_SMHM(z, mhalo):  
  '''
  Returns the median stellar mass M* for halos of mass Mh. 
  Because of scatter in the SMHM relation, the inverse of 
  this function does not give the average halo mass for a given stellar mass.

  Power law slope of -alpha for mhalo << M1. 
  Sub power-law of index gamma for mhalo >> M1.

  characterisitc stellar mass to halo mass is epsilon
  at the characteristic halo mass M1.
    
  Maximum errors are 0.025 dex, as opposed to 0.1 dex of DPL.

  Returns np.log10(<M* (mhalo)>).
  '''
  
  epsilon, alpha, delta, gamma, Sigma, M1, c = _params(z)

  x  = np.log10(mhalo / M1)
  SM = np.log10(epsilon * M1) + _f(x, alpha, delta, gamma) - _f(0.0, alpha, delta, gamma)

  return  10. ** SM, Sigma, c

## And vectorize.  
med_SMHM = np.vectorize(_med_SMHM)


if __name__ == '__main__':
    print('\n\nWelcome to the Behroozi calculator.\n\n')
    
    ## zs = np.array([0.1, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])  ## [Msun].
    zs = np.array([3.0, 4.0, 5.0])
    Ms = np.logspace(10.0, 15.0, num=50)                             ## [Msun].
    
    latexify(fig_width=None, fig_height=None, columns=1, equal=True, fontsize=10, ggplot=True, usetex=True)
    
    '''
    for z in zs:
        SMs, Errs, c = med_SMHM(z, Ms) 
        epsilon, alpha, delta, gamma, Sigma, M1, c = _params(z)
        
        pl.loglog(Ms, SMs, label=r'$z = $' + '%.1lf' % z)

        Mc = _Mc(alpha, delta, gamma, M1)

        ## pl.axvline(Mc, ymin=0., ymax=1., c='k', alpha=0.4)
    '''
    
    ##  'Implicit' [h^-1 M*] conversion. 
    for drop in ['u','g', 'r']:
      SMs, Errs, c = ishikawa(drop, Ms)
      pl.loglog(Ms * cparams['h_100'], SMs, label=r'$%s$-dropout' % drop)
    

    pl.xlabel(r'$M_{\rm{halo}} \ [h^{-1} M_{\odot}]$')
    ## pl.ylabel(r'$M_{*} \ [M_{\odot}]$')
    pl.ylabel(r'$M_{*} \ [M_{\odot}]$')   
    
    pl.xlim(1.e10, 1.e14)
    pl.ylim(1.e7,  1.e12)

    pl.legend(loc=4, ncol=1)
    
    plt.tight_layout()

    pl.savefig('plots/ishikawa.pdf')

    print('\n\nDone.\n\n')
