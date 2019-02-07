import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt
import  astropy.units      as      u 
import  astropy.constants  as      const
 
from    bz                 import  get_dropoutbz
from    prep_camb          import  CAMB
from    pmh                import  Pmm, get_PkInterps
from    params             import  get_params
from    cosmo              import  cosmo
from    scipy.integrate    import  nquad
from    growth_rate        import  growth_rate 
from    utils              import  latexify
from    goldrush.specs     import  samplestats  as goldrush_stats
from    Malkan.specs       import  samplestats  as malkan_stats
from    reddy              import  samplestats  as reddy_stats
from    schechter.nbar     import  comovdensity
from    scipy.interpolate  import  interp1d
from    completeness       import  get_dropoutpz


params = get_params()

def linb(z):
    return (1. + z)

def const_nz(z, ngal=1.e-4):
    return  ngal  ## [(h^{-1} Mpc)^-3]

def vipers_nz(z, A=3.103, z0=0.191, alpha=8.603, beta=1.448, ngal=5e3):
    ##  ngal is the expected number of galaxies per sq. deg. 
    ##  Eqn. (2) of https://arxiv.org/pdf/1303.2622.pdf
    result  =  A * ((z / z0) ** alpha) * np.exp(-(z / z0) ** beta)

    ##  Normalisation.
    dz      =                    0.01
    zz      =  np.arange(0.0, 2.0, dz)
    norm    =  dz * np.sum(A * ((zz / z0) ** alpha) * np.exp(-(zz / z0) ** beta))

    return  result / norm 

def vol_integrand(z, fsky=0.5, fkp_weighted=False, nbar=1.e-3, P0=5.e3):
    ##  Purely volume integral as a sanity check.
    ##  dV / dz [(h^{-1} Mpc)^3];  Differential comoving volume per redshift per steradian.                                                                                                                              
    ## 
    ##  Note:  cosmo.differential_comoving_volume(z) = (const.c.to('km/s') / cosmo.H(z)) * cosmo.comoving_distance(z) ** 2. = dV/dz [d\Omega] = chi^2 dChi/dz [d\Omega]. 
    dVdz = fsky * 4. * np.pi * cosmo.differential_comoving_volume(z).value * params['h_100'] ** 3.

    if fkp_weighted:
      ##  FKP volume weighting.                                                                                                                                                                                           
      nP  = nbar * P0
      fkp = nP / (1. + nP)

      return  fkp * fkp * dVdz

    else:
      return  dVdz

def Kaiser(Pk_interps, beta, z, mu, ks):
    return  (1. + beta * mu * mu)**2. * linb(z) * linb(z) * Pmm(Pk_interps, ks, z)

def nP(z, mu, k, Pk_interps, fsky, nz):
    a      = 1. / (1. + z)
    beta   = growth_rate(a) / linb(z)

    nP     = nz(z) * Kaiser(Pk_interps, beta, z, mu, k)

    return  nP / (1. + nP)

def integrand(z, mu, y, Pk_interps, fsky, nz):
    k      = np.exp(y)
    
    ##  dV / dz [(h^{-1} Mpc)^3];  Differential comoving volume per redshift per steradian.
    ##  Note:  cosmo.differential_comoving_volume(z) = (const.c.to('km/s') / cosmo.H(z)) * cosmo.comoving_distance(z) ** 2. = dV/dz [d\Omega] = chi^2 dChi/dz [d\Omega].  
    dVdz   = cosmo.differential_comoving_volume(z).value * params['h_100'] ** 3.

    a      = 1. / (1. + z) 
    beta   = growth_rate(a) / linb(z)

    ##  FKP volume weighting.
    nP     = nz(z) * Kaiser(Pk_interps, beta, z, mu, k)
    fkp    = nP / (1. + nP)
    
    ##  Effecive (S / N).
    result = fsky * np.exp(3. * y) * dVdz * fkp * fkp / np.pi

    return  result

def plot_vipers(ngal=5.e3):
    ##  Test VIPERS N(z)
    zs   = np.arange(0.0, 1.3, 0.01)
    Ns   = vipers_nz(zs, A=3.103, z0=0.191, alpha=8.603, beta=1.448, ngal=5.e3)

    pl.plot(zs, Ns, label=ngal)
    pl.xlabel(r'$z$')
    pl.ylabel(r'$dN/dz$')
    pl.legend()
    plt.tight_layout()
    pl.show()

def check_vol(fsky=0.5, fkp_weighted=False, nbar=1.e-3, P0=5.e3):
    ##  Simple volume integral as a sanity check.
    zmin    = 0.6
    zmax    = 1.2

    zranges = [[zmin, zmax]]
    args    = (fsky, fkp_weighted, nbar, P0)

    ##  vol_integrand(z, fsky=0.5, fkp_weighted=False, nbar=1.e-3, P0=5.e3)
    result  = nquad(vol_integrand, zranges, args=args, full_output=True)

    args    = (fsky, True, nbar, P0)
    fkp_wt  = nquad(vol_integrand, zranges, args=args, full_output=True)

    print('Vol: %.4lf [(h^-1 Gpc)^3], FKP Vol:  %.4lf [(h^-1 Gpc)^3], compared to Astropy: %.4lf [(h^-1 Gpc)^3]' % (result[0] / 1.e9, fkp_wt[0] / 1.e9,\
                                                                                                                    fsky * (cosmo.comoving_volume(zmax).value - cosmo.comoving_volume(zmin).value) * params['h_100'] ** 3. / 1.e9))  

def check_nP(Pk_interps, fsky=0.5, nz=const_nz):
    zs = np.arange(1.5, 6.0, 0.5)
    ks = np.logspace(-2,  0,  20)
  
    for z in zs:
      pl.semilogx(ks, nP(z, 0.0, ks, Pk_interps, fsky=fsky, nz=nz), label=str(z))

    pl.xlabel(r'$k$')
    pl.ylabel(r'$nP/(1 + nP)$')

    pl.legend()
    pl.show()


if __name__ == '__main__':
    print('\n\nWelcome to the RSD S/N calculator.')
    
    compute = True

    if compute:
        print('Loading CAMB module.')

        cambx       =  CAMB()
        Pk_interps  =  get_PkInterps(cambx)

        fsky        =  14000. / 41253.

        zmin        =  1.2
        zmax        =  1.6

        kmaxs       =  np.arange(0.1, 0.3, 0.05)
        results     =  []
    
        ##  plot_vipers(ngal=5.e3)
        ##  check_vol()

        band        =               'r'
        stats       =  goldrush_stats()

        alpha       =  stats[band]['schechter']['alpha']
        Mstar       =  stats[band]['schechter']['M_star']
        phi_star    =  stats[band]['schechter']['phi_star']

        
        drop_nz     =  lambda z: const_nz(z, ngal = 1.e-2)  ##  [(h^-1 Mpc)^-3]. 

        ##  zee, pzee =  get_dropoutpz(drop='g')
        ##  pz        =  interp1d(zee, pzee, kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)

        ##  drop_nz   =  lambda z: pz(z) * (10. ** comovdensity(z, phi_star, Mstar, alpha, type='app', mlim=25.0, band=band, printit=False))  ##  [(h_100/Mpc)^3]
        
        ## for z in np.arange(zmin, zmax, 0.01):
        ##   pl.plot(z, drop_nz(z), 'r^', markersize=2)
        '''
        zs          =  np.arange(zmin, zmax, 0.01)
        drop_nz     =  np.array([drop_nz(z) for z in zs])
    
        drop_nz     =  interp1d(zs, drop_nz, kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)
        
        pl.plot(np.arange(0., 10., 0.01), drop_nz(np.arange(0., 10., 0.01)), 'k-')

        pl.xlabel(r'$z$')
        pl.ylabel(r'$\bar n(z)$')

        plt.tight_layout()
        pl.show()
        '''

        ##  check_nP(Pk_interps, fsky=0.5, nz=drop_nz)

        for kmax in kmaxs:
            ##  y = np.log(k)
            ymin, ymax  =  np.log(5.e-2), np.log(kmax)

            ##  [[zmin, zmax], [mu min, mu max], [ymin, max]]
            ##  Note: integral symmetric in mu -> 2 * \int [0., 1.]
            ranges      = [[zmin, zmax], [0., 1.], [ymin, ymax]]
            args        = (Pk_interps, fsky, drop_nz)

            print('\nSolving for integral with kmax = %.3lf.' % kmax)
    
            result      = nquad(integrand, ranges, args=args, full_output=True)
            results.append(result[0])

            print('Solution:  %.6le' % result[0])

        results         = np.array(results)
        results         = np.sqrt(results)

        np.savetxt('dat/rsd.txt', np.c_[kmaxs, results], fmt='%.6lf \t %.6le')

    else:
        kmaxs, results  = np.loadtxt('dat/rsd.txt', unpack=True)
        
    ##  And plot ...
    latexify(columns=1, equal=True, fontsize=10, ratio=1., ggplot=True, usetex=True)

    pl.plot(kmaxs, np.log10(results))

    pl.xlabel(r'$k_{\rm{max}} \ [(h \ \rm{Mpc}^{-1})]$')
    pl.ylabel(r'$\rm{log}_{10}(S/N)$')

    pl.savefig('plots/rsd.pdf', bbox_inches='tight')
        
    print('\n\nDone.\n\n')
