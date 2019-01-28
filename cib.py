import  numpy              as      np 
import  pylab              as      pl
import  astropy.constants  as      const

from    cosmo              import  cosmo
from    utils              import  comoving_distance


def fv(nu, beta = 2., T=34., nup=4955, alpha=2.0):
    '''
    eqn. (13) of https://arxiv.org/pdf/1705.02332.pdf
    
    Input:
        nu in GHz.
    '''

    nu               = np.asarray(nu)
    exponent         = (const.h.value * nu) / (const.k_B.value * T)
    result           = nu ** (3. + beta) / (np.exp(exponent) - 1.0)

    exponent         = (const.h.value * nup) / (const.k_B.value * T)
    result[nu > nup] = nup ** (3. + beta) / (np.exp(exponent) - 1.0) / (nu[nu > nup] / nup) ** alpha

    return  result

def WCIB(z, nu = 353., zc=2., sigmaz=2.):
    '''
    Input:
        nu in GHz; [545, 857, 217, 353] GHz.
    '''

    result  = fv(nu * (1. + z), beta = 2., T=34., nup=4955, alpha=0.0)
    result *= np.exp(-0.5*(z - zc)**2. / sigmaz**2.)
    result *= comoving_distance(z) ** 2.
    result /= (cosmo.H(z).value * (1. + z)**2.)
    result /= result.max()

    return  result

def nWCIB(z, nu = 353., zc=2., sigmaz=2.):
    ##  Needs normalised.
    dz      = 0.01
    zs      = np.arange(0.0, 10., dz)

    Ws      = WCIB(zs, nu, zc, sigmaz)
    norm    = dz * np.sum(Ws)

    return  WCIB(z, nu, zc, sigmaz) / norm

def cib_effshot():
    ##  353GHz CIB
    return  2.1e-13


if __name__ == '__main__':
    print('\n\nWelcome to CIB.\n\n')
    
    zs  = np.arange(0.0, 10.0, 0.01)
    Ws  = WCIB(zs) 

    pl.plot(zs, Ws)

    pl.savefig('plots/WCIB.pdf')

    print('\n\nDone.\n\n')
