import  hmf
import  numpy              as      np

from    growth_rate        import  growth_factor 
from    pmh                import  Pmm, get_PkInterps, linz_bz
from    prep_camb          import  CAMB, Clxy
from    scipy.integrate    import  simps
from    params             import  get_params
from    scipy.interpolate  import  interp1d
from    tinker             import  tinker_fsigma
from    hmf                import  MassFunction


params         = get_params()
params['om_m'] = 0.30851

rho_b          = 2.78e11 * params['om_m'] * params['h_100'] ** 2. / params['h_100'] ** 3.   ## [Msol * (Mpc/h)^3] 


## Table 3 of https://arxiv.org/pdf/1005.2239.pdf
def jenkins(sigma):
    ## Valid for 0 < z < 5.
    arg      = - np.log(sigma)

    exponent = - np.abs(arg + 0.61) ** 3.8
    result   = 0.315 * np.exp(exponent)
    
    print('Jenkins mass fn. limited to -1.2 < ln |1/sigma| <= 1.05 and 0 < z < 5.')

    result[arg < -1.20] = np.nan
    result[arg >  1.05] = np.nan

    return  result

def sheth_tormen(sigma):
    deltac  = 1.686    
    
    fsigma  = 0.3222 * np.sqrt( 2. * 0.75 / np.pi) * np.exp( - 0.75 * deltac ** 2. / 2. / sigma ** 2.) 
    fsigma *= (1. + (sigma**2. / 0.75 / deltac ** 2.) ** 0.3) * deltac / sigma
    
    return  fsigma
    
def W(x):
    return 3. * (np.sin(x) - x * np.cos(x)) / x**3.

def sigma(Pk_interps, Ms, z):
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
      integrand = ks * ks * Ps * Ws * Ws 

      sig2      = simps(integrand, dx = ks[1] - ks[0])

      sig2     *= growth_factor(z) ** 2.
      sig2     /= (2. * np.pi ** 2.)
      
      result.append(np.sqrt(sig2))

    return  np.array(result)

def tinker_bias(nu, Delta=200):
    ## Tinker bias fn. (https://arxiv.org/pdf/1001.3162.pdf)
    y = np.log10(Delta)

    A = 1.0 + 0.24 * y * np.exp(- (4/y)**4.)
    a = 0.44 * y - 0.88
    B = 0.183
    b = 1.5
    C = 0.019 + 0.107 * y + 0.19 * np.exp(-(4./y)**4.)
    c = 2.4

    return  1. - A * nu ** a /(nu ** a + 1.686 ** a) + B * nu ** b + C * nu **c


if __name__ == "__main__":
    import  matplotlib         as      mpl
    import  matplotlib.pyplot  as      plt
    import  pylab              as      pl
    import  matplotlib.patches as      patches

    from    prep_camb          import  CAMB, Clxy
    ## from    nbar               import  comovdensity
    ## from    specs              import  samplestats
    from    pylab              import  rcParams


    ## plt.style.use('ggplot')

    mpl.rc('text', usetex = True)

    rcParams['figure.figsize'] = (3.5, 3.5)

    print("\n\nWelcome to massfn.\n\n")
    
    z                   =  4.0
    Ms                  =  1.e9 * 10.** np.arange(0.0, 5.8, 0.05)
    '''
    ##  Prepare pycamb module; linear, non-linear matter P(k) and Cls.                                                                                      
    cambx               =  CAMB()
    Pk_interps          =  get_PkInterps(cambx)

    for lsigma, lfsigma, label in zip([sigma, sigma, sigma], [sheth_tormen, jenkins, tinker_fsigma], ['ST', 'Jenkins', 'Tinker']):
        sigmas          =  lsigma(Pk_interps, Ms, z)
        fsigmas         =  lfsigma(sigmas)

        ## nus          =  1.686 / sigmas
        ## tinker_bs    =  tinker_bias(nus, Delta=200)

        dndsig          =  - rho_b * fsigmas / Ms / sigmas
    
        ns              = []
        ## Meff         = []

        for M in Ms[:-1]:
            ns.append(simps(-dndsig[Ms > M], dx = np.abs(sigmas[1] - sigmas[0])))
            ## Meff.append( np.abs( simps(-dndsig[Ms > M] * Ms[Ms > M][::-1], dx = np.abs(sigmas[1] - sigmas[0])) ) )
        
        ns       = np.array(ns)
        ## Meff  = np.array(Meff) / ns
        
        pl.loglog(Ms[:-1], ns, label=label)
    '''
    ## Load Martin's results.
    data = np.loadtxt('../dat/massfn_st_0.2000.dat')
    pl.loglog(data[:,0], data[:,1], label='Martin-ST')

    data = np.loadtxt('../dat/massfn_tin_0.2000.dat')
    pl.loglog(data[:,0], data[:,1], label='Martin-Tinker')
    
    ##  Murray et al. hmf python code. 
    mf   = MassFunction()

    for label, fitting, color in zip(['ST', 'Tinker'], [hmf.fitting_functions.ST, hmf.fitting_functions.Tinker08], ['c', 'r']):
        mf.update(z = 4.0, Mmin = 8., Mmax = np.log10(5.e14), cosmo_params = {'Om0': 0.30841}, hmf_model = fitting, delta_h = 200.)

        print  mf.parameter_values
        print  mf.parameter_info()

        pl.loglog(mf.m, mf.ngtm, label='py-hmf: ' + label, c=color, alpha=0.6)


    ## Axis labels.
    ax   = pl.gca()

    ax.set_xlabel(r'$M_{\rm{min}} \ [M_\odot / h]$')
    ax.set_ylabel(r'$\bar n \ [(h^{-1} \ \rm{Mpc})^{-3}]$')

    pl.ylim(1.e-10, 1.e2)

    pl.legend(loc=3)

    pl.savefig('../plots/massfn.pdf', bbox_inches='tight')
    
    print("\n\nDone.\n\n")
