import  os 
import  numpy              as      np
import  pylab              as      pl
import  matplotlib         as      mpl
import  matplotlib.pyplot  as      plt 

from    nbar               import  mlimitedM, dVols
from    params             import  get_params
from    qso_kcorr          import  get_qsokcorr
from    cosmo              import  cosmo


params = get_params()

def  get_Mgs(zee):
    ##  Table 4 of https://arxiv.org/pdf/1509.05607.pdf                                                                                                   
    zp          =     2.2
    Mgs_zp      =  -26.71

    ##  Parameters split across zpivot.                                                                                                                    
    if (0.68 <= zee) & (zee <= 2.2):
      alpha       =  -4.31
      beta        =  -1.54

      k1          =  -0.08
      k2          =  -0.40

    elif (2.2 < zee) & (zee <= 4.0):
      alpha       =  -3.04
      beta        =  -1.38

      k1          =  -0.25
      k2          =  -0.05

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
    if (0.68 <= zee) & (zee <= 2.2):
      alpha       =  -4.31
      beta        =  -1.54

    elif (2.2 < zee) & (zee < 4.0):
      alpha       =  -3.04
      beta        =  -1.38

    else:
        raise ValueError("\n\nQSO luminosity fn. of Palanque-Delabrouille is not defined at this redshift: %.3lf" % zee)

    Phi         = Phi_star
    denom       = 10. ** (0.4 * (1. + alpha) * (Mg - Mgs)) + 10. ** (0.4 * (1. + beta) * (Mg - Mgs)) 

    Phi        /=    denom

    return  Phi

def kcorr(zee):
    if zee > 4.0:
        raise ValueError('NPD k-correction (McGreer) not defined for z<0.6 or z>4.')

    zs, Ks, Ksp = get_qsokcorr(plotit=False)

    return  Ksp(zee)

def gmag(Mg, zee, restM=False, printit=False):
    ##  Observationally, this is the galactic extinction corrected g band mag. 
    ##  i.e. gdered of NPD++.
    ##  Mg is the z=2 absolute mag. in the g-band. 
    ##  eqn. (4) of https://arxiv.org/pdf/1509.05607.pdf

    ##  g-band k-correction.
    k = kcorr(zee)

    if not restM:
        ##  If not rest Absolute Mag, then assumed z=2 as per NPD.
        k -= kcorr(2.0)

    if printit:
      print(zee)
      print(Mg)
      print(kcorr(zee), kcorr(2.0))
      print(cosmo.distmod(zee))

    return  Mg + cosmo.distmod(zee).value + k

def get_ns(Ms, zee=2.5):
    Phis  = get_Phi(Ms, zee)

    dM    =  Ms[1] - Ms[0]

    ns    =  np.cumsum(Phis) * dM
    ns   /=  params['h_100'] ** 3.

    return  ns


if __name__ == '__main__':
    print('\n\nWelcome to a qso luminosity fn. calculator.\n\n')

    dM, dz = 0.1, 0.4
    zs     = np.arange( 0.68,   3.2, dz)
    Ms     = np.arange(-32.,   -15., dM)

    result = []

    for zee in zs:        
      Phis  = get_Phi(Ms, zee)
      ns    = get_ns(Ms, zee=zee)

      pl.semilogy(Ms, Phis, label=r'$z=$' + '%.1lf' % zee)
    
    ## Fig 10. of https://arxiv.org/pdf/1509.05607.pdf  
    pl.xlabel(r'$M_g(z=2)$')
    pl.ylabel(r'$\Phi(M_g, z) = dN/dM/dV \ [\rm{Mpc}^{-3} \rm{\ per \ mag}]$')
    
    pl.legend(ncol=2)
    pl.show()

    ##  Number counts. 
    pl.clf()

    for zee in zs:
      ## Mlo to Mhi -> faintest object included in bin is M+dM.
      gs = gmag(Ms + dM, zee, restM=False, printit=False)
      ns = get_ns(Ms, zee=zee)

      pl.semilogy(gs, ns, label=r'$z=$' + '%.1lf' % zee)

    pl.xlim(20., 26.)
    pl.ylim(1.e-6, 1.e-2)

    pl.xlabel(r'$g_{AB}$')
    pl.ylabel(r'$n$ [(h/Mpc)$^3$]')

    pl.legend(ncol=2)
    pl.show()

    print('\n\nDone.\n\n')
