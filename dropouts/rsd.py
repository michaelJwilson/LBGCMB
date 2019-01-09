import  numpy            as      np
import  pylab            as      pl

from    bz               import  get_dropoutbz
from    prep_camb        import  CAMB
from    pmh              import  Pmm, get_PkInterps
from    params           import  get_params
from    cosmo            import  cosmo
from    scipy.integrate  import  nquad
from    growth_rate      import  growth_rate 
from    utils            import  latexify


def linb(z):
    return (1. + z)

def Gauss_nz(z, z0=2., sigz=0.5):
    return  np.exp(-(z-z0)**2. / 2. / sigz**2.) / np.sqrt(2. * np.pi * sigz**2.)

def Kaiser(Pk_interps, beta, z, mu, ks):
    return  (1. + beta * mu * mu)**2. * linb(z) * linb(z) * Pmm(Pk_interps, ks, z) 

def integrand(z, mu, y, Pk_interps, fsky):
    k    = np.exp(y)
    
    ##  dV/dz [(h^{-1} Mpc)^3]
    dVdz = fsky * 4. * np.pi * cosmo.differential_comoving_volume(z).value * params['h_100'] ** 3.

    a    = 1. / (1. + z) 
    beta = growth_rate(a) / linb(z)

    ## FKP volume weighting
    nP   = Gauss_nz(z) * Kaiser(Pk_interps, beta, z, mu, k)
    fkp  = nP / (1. + nP)
    
    ##  Effecive S / N.
    return  np.exp(3. * y) * dVdz * fkp * fkp / (2. * np.pi) 


if __name__ == '__main__':
    print('\n\nWelcome to the RSD S/N calculator.')

    compute     =  True
    params      =  get_params()

    if compute:
        print('Loading CAMB module.')

        cambx       =  CAMB()
        Pk_interps  =  get_PkInterps(cambx)

        fsky        =  0.5
        zmin        =  2.0
        zmax        =  3.0

        kmaxs       =  np.arange(0.2, 0.45, 0.05)
        results     =  []
    
        for kmax in kmaxs:
            ##  y = np.log(k)
            ymin, ymax  =  np.log(1.e-3), np.log(kmax)

            ## [[zmin, zmax], [mu min, mu max], [ymin, max]]
            ## Note: integral symmetric in mu -> 2 * \int [0., 1.]
            ranges      = [[zmin, zmax], [0., 1.], [ymin, ymax]]
            args        = (Pk_interps, fsky)

            print('Solving for integral with kmax = %.4le.' % kmax)
    
            result      = nquad(integrand, ranges, args=args, full_output=True)
            results.append(result[0])

            print('Solution for integral:', result)

        results = np.array(results)
        results = np.sqrt(results)

    else:
        kmaxs, results = np.loadtxt('dat/rsd.txt', unpack=True)
        

    latexify(columns=2, equal=False, fontsize=10, ratio=None, ggplot=True, usetex=True)

    pl.plot(kmaxs, np.log10(results))

    pl.xlabel(r'$k_{\rm{max}} \ [(h \ \rm{Mpc}^{-1})]$')
    pl.ylabel(r'$\rm{log}_{10}(S/N)$')

    pl.savefig('plots/rsd.pdf', bbox_inches='tight')
    
    np.savetxt('dat/rsd.txt', np.c_[kmaxs, results], fmt='%.6lf \t %.6le')

    print('\n\nDone.\n\n')
