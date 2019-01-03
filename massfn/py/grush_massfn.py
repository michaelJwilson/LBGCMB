import  numpy              as      np

from    growth_rate        import  growth_factor 
from    pmh                import  Pmm, get_PkInterps, linz_bz
from    prep_camb          import  CAMB, Clxy
from    scipy.integrate    import  simps
from    params             import  get_params
from    scipy.interpolate  import  interp1d
from    tinker             import  tinker_sigma, tinker_fsigma


params         = get_params()
params['om_m'] = 0.30851

rho_b          = 2.78e11 * params['om_m'] * params['h_100'] ** 2. / params['h_100'] ** 3.   ## [Msol * (Mpc/h)^3] 


## Table 3 of https://arxiv.org/pdf/1005.2239.pdf
def jenkins(sigma):
    ## Valid for 0 < z < 5.
    arg      = np.log(1. / sigma)

    exponent = - np.abs(arg + 0.61) ** 3.8
    result   = 0.315 * np.exp(exponent)
    
    print('Jenkins mass fn. limited to -1.2 < ln |1/sigma| <= 1.05 and 0 < z < 5.')

    result[arg < -1.20] = np.nan
    result[arg >  1.05] = np.nan

    return  result

def sheth_tormen(sigma):
    deltac = 1.686    
    fsigma = 0.3222 * np.sqrt( 2. * 0.75 / np.pi) * np.exp( - 0.75 * deltac ** 2. / 2. / sigma ** 2.) * (1. + (sigma**2. / 0.75 / deltac ** 2.) ** 0.3) * deltac / sigma
    
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
    Ms                  =  1.e11 * 10.** np.arange(0.0, 3.8, 0.05)
    
    ##  Prepare pycamb module; linear, non-linear matter P(k) and Cls.                                                                                      
    cambx               =  CAMB()
    Pk_interps          =  get_PkInterps(cambx)

    for lsigma, lfsigma, label in zip([sigma, sigma, tinker_sigma], [sheth_tormen, jenkins, tinker_fsigma], ['ST', 'Jenkins', 'Tinker']):
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
    band     = 'g'
    stats    = samplestats()

    phi_star = stats[band]['schechter']['phi_star'] 
    M_star   = stats[band]['schechter']['M_star']
    alpha    = stats[band]['schechter']['alpha']
    
    ## Second y-axis with app. mag. limits for GoldRush. 
    mags     = np.arange(10.0, 35.0, 0.5)
    nbars    = []
    
    for mag in mags:
        nbars.append( 10. ** comovdensity(z, phi_star, M_star, alpha, type='app', mlim=mag, printit=False) )

    interp_nbar = interp1d(nbars, mags, bounds_error=True)
    mags_ax     = interp_nbar(ns)
    '''
    ax          = pl.gca()
    '''
    ## Add second mag_axis.
    ax2         = ax.twinx()

    ax2.semilogx(Ms[:-1], mags_ax, alpha=0.0)
    ax2.set_ylabel(r'$i_{{\rm{AB}}}$')
    
    ## Add b_eff axis on top.
    ax3         = ax.twiny()

    ax3.semilogy(tinker_bs[:-1], ns, alpha=0.0)
    ax3.set_xlabel(r'$b_{\rm{min}}$')
    '''
    ## and original
    ax.set_xlabel(r'$M_{\rm{min}} \ [M_\odot]$')
    ax.set_ylabel(r'$\bar n \ [(h^{-1} \ \rm{Mpc})^{-3}]$')
    '''
    ## Add arrows
    style = 'Simple, tail_width=0.5, head_width=4, head_length=8'
    kw    = dict(arrowstyle=style, color='k', alpha=0.3)

    ## Schechter fn.                                                                                                                                       
    arrow = patches.FancyArrowPatch((ax.get_xlim()[1], 1.e-3), (ax.get_xlim()[0], 1.e-3), connectionstyle = 'arc3, rad=-.3', **kw)

    ax.add_patch(arrow)
    ax.annotate('Schechter fn.', xy=(3.e14, 4.e-5), xycoords='data', va='bottom', ha='left') ## bbox=dict(fc='w') 

    ## Jenkins mass fn.
    arrow = patches.FancyArrowPatch((ax.get_xlim()[0], 1.e-3), (1.e13,   ax.get_ylim()[0]), connectionstyle = 'arc3, rad=-.3', **kw)

    ax.add_patch(arrow)
    ax.annotate('Jenkins\n(2001)', xy=(9.e11, 5.e-6), xycoords='data', va='bottom', ha='left') ## bbox=dict(fc='w')
 
    ##  Tinker bias.
    arrow = patches.FancyArrowPatch((1.e13, ax.get_ylim()[0]), (1.e14,   ax.get_ylim()[1]), connectionstyle = 'arc3, rad=-.3', **kw)

    ax.add_patch(arrow)
    ax.annotate('Tinker\n(2010)', xy=(1.4e12, 1.e-7), xycoords='data', va='bottom', ha='left') ## bbox=dict(fc='w')
    
    ## Harikane 17.
    kw    = dict(arrowstyle=style, color='g')
    arrow = patches.FancyArrowPatch((ax.get_xlim()[1], 1.e-3), (1.5e13,  ax.get_ylim()[1]), connectionstyle = 'arc3, rad=-.3', **kw)

    ax.add_patch(arrow)
    ax.annotate('Harikane (2017)', xy=(2.5e14, 1.3e-3), xycoords='data', va='bottom', ha='left') ## bbox=dict(fc='w')

    ## Straight line connection.
    ax.axhline(y = 1.e-3, xmin = ax.get_xlim()[0], xmax = ax.get_xlim()[1], c='k')
    '''
    pl.legend(loc=1)

    pl.savefig('../plots/abmatch_linb.pdf', bbox_inches='tight')
    
    print("\n\nDone.\n\n")
