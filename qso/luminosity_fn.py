import  numpy              as      np
import  pylab              as      pl
import  matplotlib         as      mpl
import  matplotlib.pyplot  as      plt 

from    nbar               import  mlimitedM, dVols
from    params             import  get_params


plt.style.use('ggplot')
mpl.rc('text', usetex = True)

params = get_params()

def  get_Mgs(zee):
    ##  Table 4 of https://arxiv.org/pdf/1509.05607.pdf                                                                                                   
    zp          =     2.2
    Mgs_zp      =  -26.71

    ##  Parameters split across zpivot.                                                                                                                    
    if (0.68 < zee) & (zee < 2.2):
      alpha       =  -4.31
      beta        =  -1.54

      k1          =  -0.08
      k2          =  -0.40

    elif (2.2 < zee) & (zee < 4.0):
      alpha       =  -4.31
      beta        =  -1.54

      k1          =  -0.08
      k2          =  -0.40

    else:
        raise ValueError("\n\nQSO luminosity fn. of Palanque-Delabrouille is not defined at this redshift: %.3lf" % zee)

    ##  See eqn. (7) of https://arxiv.org/pdf/1509.05607.pdf                                                                                               
    Mgs = Mgs_zp - 2.5 * ( k1 * (zee - zp) + k2 * (zee - zp) ** 2. )

    return  Mgs

def get_Phi(Mg, zee):
    zp          =     2.2
    log_Phistar =  -6.010
    
    Phi_star    = 10. ** log_Phistar

    Mg          = np.asarray(Mg)
    Mgs         = get_Mgs(zee)

    ##  Parameters split across zpivot.                                                                                                                    
    if (0.68 < zee) & (zee < 2.2):
      alpha       =  -4.31
      beta        =  -1.54

    elif (2.2 < zee) & (zee < 4.0):
      alpha       =  -4.31
      beta        =  -1.54

    else:
        raise ValueError("\n\nQSO luminosity fn. of Palanque-Delabrouille is not defined at this redshift: %.3lf" % zee)

    Phi         = Phi_star
    denom       = 10. ** (0.4 * (1. + alpha) * (Mg - Mgs)) + 10. ** (0.4 * (1. + beta) * (Mg - Mgs)) 

    Phi        /=    denom

    return  Phi

def kcorr(zee):
    alphav = -0.5

    return  -2.5 * (1. + alphav) * np.log10(1. + zee)

def get_ns(Ms, zee=2.5):
    Phis  =  get_Phi(Ms, zee)

    dM    =  Ms[1] - Ms[0]

    ns    =  np.cumsum(Phis * dM)
    ns   /=  params['h_100'] ** 3.

    return  ns


if __name__ == '__main__':
    print('\n\nWelcome to a qso luminosity fn. calculator.\n\n')

    zee = 2.5

    dm  =  0.1

    ms  = np.arange(22.5, 25.5, dm)
    Ms  = mlimitedM(zee, ms, M_standard=None, kcorr=True)

    ns  =  get_ns(Ms)

 
    pl.semilogy(ms, ns, label=r'$2.2 \ \leq \ z \ \leq \ 2.6$')

    pl.xlabel(r'$m_g$')
    pl.ylabel(r'$n(m < m_g) \quad [(h^{-1} \ \rm{Mpc})^{-3}]$')

    pl.xlim(-29.,   -15.)
    pl.ylim(5.e-8, 1.e-3)
    
    pl.legend(loc=2)

    pl.savefig('plots/qso_lumfn.pdf')

    '''
    pl.clf()

    ## Plot simple k-correction. 
    zs  = np.arange(0.0, 4.5, 0.01)
    ks  = kcorr(zs)

    pl.plot(zs, ks, 'r-')

    dat = np.loadtxt('kcorr/table4.dat')
    pl.plot(dat[:,0], dat[:,1], 'k-')

    pl.xlim( 0.5,  3.0)
    pl.ylim( 0.2, -1.4)

    pl.xlabel(r'$z$')
    
    pl.savefig('plots/kcorr.pdf')
    '''
    print('\n\nDone.\n\n')
