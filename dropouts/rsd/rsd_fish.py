import  os
import  sys
import  vegas
import  pickle
import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt
import  astropy.units      as      u 
import  astropy.constants  as      const
 
from    prep_camb          import  CAMB
from    pmh                import  Pmm, get_PkInterps
from    params             import  get_params
from    cosmo              import  cosmo
from    scipy.integrate    import  nquad
from    growth_rate        import  growth_rate 
from    utils              import  latexify
from    scipy.interpolate  import  interp1d
from    numpy.linalg       import  inv
from    euclid             import  euclid_bz, euclid_nz, euclid_area
from    schechter.gen_pz   import  peakz           as  _peakz
from    schechter.get_shot import  get_shot
from    schechter.get_pz   import  get_pz
from    nbar               import  comovdensity
from    get_bz             import  bz_callmodel
from    get_schechters     import  get_schechters
from    get_wschechters    import  get_wschechter
from    reddy.specs        import  samplestats       as  reddy_stats
from    goldrush.specs     import  samplestats       as  grush_stats
from    Malkan.specs       import  samplestats       as  malkan_stats
from    reddy.pz           import  get_pz            as  reddy_getpz
from    reddy.completeness import  get_completeness  as  reddy_getcompleteness
from    Malkan             import  completeness      as  malkan_completeness
from    goldrush           import  completeness      as  grush_completeness


cparams = get_params()

def Kaiser(Pk_interps, beta, z, mu, ks, bz, type='linear', sig=0.0):  
    return  (1. + beta * mu * mu) ** 2. * bz(z) * bz(z) * Pmm(Pk_interps, ks, z, type=type) * np.exp(-ks * ks * mu * mu * sig * sig)

def fish_weight(b, f, mu, k, coeff):
    if coeff   == 'b':
        return   2. / (b + f * mu ** 2.)

    elif coeff == 'f':
        return   2. * mu * mu / (b + f * mu ** 2.)

    elif coeff == 's':
        return  -k * k * mu * mu 

    else:
        raise  ValueError('Unacceptable coeff: %s' % coeff)

def integrand(z, mu, y, Pk_interps, fsky, nz, bz, fish_coeff=None, type='linear', deltav=0.0):
    k      = np.exp(y)
    
    ##  dV / dz [(h^{-1} Mpc)^3];  Differential comoving volume per redshift per steradian.
    ##  Note:  cosmo.differential_comoving_volume(z) = (const.c.to('km/s') / cosmo.H(z)) * cosmo.comoving_distance(z) ** 2. 
    ##                                               = dV/dz [d\Omega] = chi^2 dChi/dz [d\Omega].  
    ## 

    dVdz   = cosmo.differential_comoving_volume(z).value * cparams['h_100'] ** 3.

    a      = 1. / (1. + z) 
    f      = growth_rate(a)

    b      = bz(z)
    beta   = f / b

    sigp   =  (1. + z) * deltav  / cosmo.efunc(z) / 100.  ##  [Mpc / h]. 

    ##  FKP volume weighting.
    nP     = nz(z) * Kaiser(Pk_interps, beta, z, mu, k, bz, type=type, sig=sigp)
    fkp    = nP / (1. + nP)
    
    ##  Effecive (S / N).
    result = np.exp(3. * y) * fsky * dVdz * fkp * fkp / np.pi

    if fish_coeff == None:
        return  result

    else:
        return  result * fish_weight(b, f, mu, k, fish_coeff[0]) * fish_weight(b, f, mu, k, fish_coeff[1])

def _vintegrand(x, args):
    ##  Vegas wrapper of integrand; input args as a list. 
    z, mu, y                                             = x[0], x[1], x[2]
    (Pk_interps, fsky, nz, bz, fish_coeff, type, deltav) = args 

    return  integrand(z, mu, y, Pk_interps, fsky, nz, bz, fish_coeff, type, deltav)


if __name__ == '__main__':
    print('\n\nWelcome to the RSD S/N calculator.')

    print('Loading CAMB module.')

    band        =     'BX'
    area        =   15000. 
        
    deltav      =     400.    ##  [km / s].
    kmax        =     0.9     ##  Maximum wavenumber tested. 

    ##
    setup       = {'BX': {'maglim': 24.0, 'decband': 'R', 'specfrac': 1.0,  'stats': reddy_stats(),  'CC':  reddy_getcompleteness()},\
                    'u': {'maglim': 24.0, 'decband': 'R', 'specfrac': 1.0,  'stats': malkan_stats(), 'CC':  malkan_completeness.get_completeness()},\
                    'g': {'maglim': 25.5, 'decband': 'i', 'specfrac': 0.1,  'stats': grush_stats(),  'CC':  grush_completeness.get_completeness('g')},\
                    'r': {'maglim': 25.5, 'decband': 'z', 'specfrac': 0.1,  'stats': grush_stats(),  'CC':  grush_completeness.get_completeness('r')}}
        
    CC          =  setup[band]['CC']
    mlim        =  setup[band]['maglim']
    stats       =  setup[band]['stats']
    specfrac    =  setup[band]['specfrac']

    ##  midz, alpha, M_star, phi_star = get_wschechter(4.0)  
    midz, alpha, M_star, phi_star     = get_schechters(stats, band) 
    
    bz          =  lambda z:  bz_callmodel(z, mlim)
    nz          =  lambda z:  specfrac * CC(z) * 10. ** comovdensity(z, phi_star, M_star, alpha, type='app', mlim=mlim, printit=False)

    pz          =  get_pz(band)
    peakz       =    _peakz(pz)

    zmin        =  peakz - 1.0
    zmax        =  peakz + 1.0

    fsky        =  area / 41253.

    dkmax       =  0.1
    kmaxs       =  np.arange(0.1, kmax + dkmax, dkmax)

    ## 
    cambx       =  CAMB()
    Pk_interps  =  get_PkInterps(cambx)

    params      =  ['b', 'f']       ##  ['b', 's', 'f']                                                                                                 
    results     =  []

    for kmax in kmaxs:
        ##  y = np.log(k)
        ymin, ymax  =  np.log(5.e-2), np.log(kmax)

        ##  [[zmin, zmax], [mu min, mu max], [ymin, max]]
        ##  Note: integral symmetric in mu -> 2 * \int [0., 1.]
        ranges      =  [[zmin, zmax], [0., 1.], [ymin, ymax]]

        integ       =  vegas.Integrator(ranges)
        
        Fisher      =  np.zeros(len(params) * len(params)).reshape(len(params), len(params))
                                    
        for i, b in enumerate(params):
            for j, f in enumerate(params):
                args    = (Pk_interps, fsky, nz, bz, b + f, 'nlinear', deltav)
                result  = integ(lambda x: _vintegrand(x, args), nitn=10, neval=1000)

                print(result.summary())

                Fisher[i, j] = result.mean
            
        ##  invert ... 
        iFisher = inv(Fisher)
            
        print(kmax, band, zmin, zmax, bz(peakz), growth_rate(1. / (1. + peakz)))
        print(kmax, band, zmin, zmax, 100. * np.sqrt(iFisher[0,0]) / bz(peakz), 100. * np.sqrt(iFisher[-1,-1]) / growth_rate(1. / (1. + peakz)))

        results.append(iFisher[-1,-1])

    results = np.array(results)
    results =  np.sqrt(results)

    print(results)

    np.savetxt('dat/rsd_%s_%.1lf_%.1lf_%.3lf.dat' % (band, kmax, fsky, deltav), np.c_[kmaxs, results], fmt='%.6lf \t %.6le')
        
    print('\n\nDone.\n\n')
