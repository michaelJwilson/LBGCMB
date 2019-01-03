import  numpy              as      np

from    scipy.integrate    import  simps
from    scipy.interpolate  import  interp1d
from    params             import  get_params
from    growth_rate        import  growth_factor


params         = get_params()
params['om_m'] = 0.30851

rho_b          = 2.78e11 * params['om_m'] * params['h_100'] ** 2. / params['h_100'] ** 3.   ## [Msol * (Mpc/h)^3]                                                                                       
def W(x):
    return 3. * (np.sin(x) - x * np.cos(x)) / x**3.

def tinker_fsigma(sigma, triangle=200):
    ## https://arxiv.org/pdf/0803.2706.pdf
    params  = {200: {'A': 0.186, 'a': 1.47, 'b': 2.57, 'c': 1.19},\
               300: {'A': 0.200, 'a': 1.52, 'b': 2.25, 'c': 1.27},\
               400: {'A': 0.212, 'a': 1.56, 'b': 2.05, 'c': 1.34},\
               600: {'A': 0.218, 'a': 1.61, 'b': 1.87, 'c': 1.45}}

    fsigma  = params[triangle]['A'] * (1. + (sigma / params[triangle]['b'] ** -params[triangle]['a']))
    fsigma *= np.exp(- params[triangle]['c'] / sigma ** 2.)

    return  fsigma
    
def tinker_sigma(Pk_interps, Ms, z):
    ## Beware variable definitions wrt massfn.py, which is based upon https://arxiv.org/pdf/1005.2239.pdf
    result = []

    ## Lift Martin's linear z=0 P(k).                                                                                                     
    data   = np.loadtxt('../dat/pklin_1.0000.txt')
    PMW    = interp1d(data[:,0], data[:,1], kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)

    for M in Ms:
      ## Implicitly a fn. of z by rho_b, Pmm and the growth factor.                                                                        
      ks        = np.arange(1.e-3, 1.e1, 1.e-3)
      
      Rms       = (3. * M / 4 / np.pi / rho_b) ** (1. / 3.)

      ## Ps     = Pmm(Pk_interps, ks, 0.0, type='linear')                                                                                   
      Ps        = PMW(ks)

      Ws        = W(ks * Rms)

      ## Defined by linear or non-linear.                                                                                            
      integrand = ks * ks * Ps * Ws

      sig       = simps(integrand, dx = ks[1] - ks[0])
      sig      *= growth_factor(z)

      result.append(sig)

    return  np.array(result)
