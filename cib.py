import  numpy              as      np 
import  pylab              as      pl
import  astropy.constants  as      const
import  matplotlib.pyplot  as      plt

from    cosmo              import  cosmo
from    utils              import  comoving_distance, latexify


def fv(nu, beta = 2., T=34., nup=4955, alpha=2.0):
    '''
    eqn. (13) of https://arxiv.org/pdf/1705.02332.pdf
    
    Input:
        nu in GHz.
    '''

    nu               =  np.asarray(nu)
    exponent         = (const.h.value * nu) / (const.k_B.value * T)
    result           =  nu ** (3. + beta) / (np.exp(exponent) - 1.0)

    exponent         = (const.h.value * nup) / (const.k_B.value * T)
    result[nu > nup] =  nup ** (3. + beta) / (np.exp(exponent) - 1.0) / (nu[nu > nup] / nup) ** alpha

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

def cib_shot():
    ##  353 GHz CIB
    return  2.1e-13


if __name__ == '__main__':
  latexify(fig_width=None, fig_height=None, columns=1, equal=True, fontsize=10)

  print('\n\nWelcome to CIB.\n\n')
    
  zs  = np.arange(0.0, 10.0, 0.01)
  Ws  = WCIB(zs) 

  pl.plot(zs, Ws)

  pl.xlabel(r'$z$')
  pl.ylabel(r'$W_{\rm{CIB}}$')

  plt.tight_layout()

  pl.savefig('plots/WCIB.pdf')
  
  print('\n\nDone.\n\n')
