import  pylab         as      pl
import  numpy         as      np

from    scipy.special import  gamma, gammaincc
from    cosmo         import  cosmo
from    params        import  get_params 
from    nbar          import  dVols
from    massfn        import  tinker_bias


params = get_params()

def Sobral_LyaSch(z):
    if (z < 2.2) | (z > 5.4):
        raise  ValueError('Sobral SC4k defined for 2.2 < z < 5.4 only.')

    ##  Table 6 of https://arxiv.org/pdf/1712.04451.pdf
    ##  redshift    alpha    log10(L*) [erg/s]    log10(Phi*) [Mpc ** -3]    rho Lya /10^40 [ergs/s/Mpc3].
    sc4kfits =  np.array([[2.2, -1.8, 42.69, -3.33, 0.48], [2.5, -1.80, 42.76, -3.23, 0.73], [3.1, -1.80, 42.69, -2.73, 1.90], [3.9, -1.80, 42.89, -3.71, 0.34], [4.7, -1.80, 43.10, -3.82, 0.48], [5.4, -1.80, 43.35, -4.18, 0.41]])
    synfits  =  np.array([[2.2, -2.0, 42.82, -3.59, 0.52], [2.5, -1.72, 42.71, -3.10, 0.74], [3.1, -1.63, 42.77, -3.06, 0.86], [3.9, -2.26, 42.93, -3.66, 1.11], [4.7, -2.35, 43.28, -4.25, 1.16], [5.4, -1.98, 43.28, -3.83, 1.11]])

    ##  Find available redshift closest to requested.
    ind      =  np.where(np.abs(synfits[:,0] - z) == np.abs(synfits[:,0] - z).min())

    alpha    =         synfits[ind,1]
    L_star   =  10. ** synfits[ind,2]                          ##  [ergs / s]
    phi_star =  10. ** synfits[ind,3] / params['h_100'] ** 3.  ##  [(Mpc / h)^-3]  

    return  alpha, L_star, phi_star

def logMhalo(logLya, z):
  ##  Eqn. (13) of https://arxiv.org/pdf/1811.00556.pdf                                                                                                                                                                                    
  ##  Note implicit dependence on redshift via L*.
  
  alpha, L_star, phi_star = Sobral_LyaSch(z)

  Lya = 10. ** logLya

  if Lya < L_star:
    Mhalo = 10. ** 12.19 * (Lya / L_star) ** 2.08  ## [Msun / h]  

  else:  
    Mhalo = 10. ** 12.19 * (Lya / L_star) ** 0.63  ## [Msun / h]  

  return  np.log10(Mhalo)

def lya_nbar(z, logLmin=None, printit=False):
    alpha, L_star, phi_star = Sobral_LyaSch(z)

    if logLmin == None:
      ##  Set lower limit on luminosity integral of L*.   
      Lmin   =  L_star

    else:
      Lmin   =  10. ** logLmin

    wmin     =  np.log10(Lmin / L_star)

    dw       =  0.1
    ws       =  np.arange(wmin, wmin + 4., dw)  ##  w = log10(L / L*)
    xs       =  np.log(10) * ws                 ##  x =    ln(L / L*) 
    ys       =  10. ** ws                       ##  y =       L / L*  

    nbar     =  phi_star * np.exp(-ys) * 10.**((1. + alpha) * ws)
    nbar     =  np.log(10.) * np.sum(nbar) * dw                 

    if printit:
      print('%.3lf  %.3lf  %.3lf  %.3le  %.3le  %.3le' % (z, synfits[ind,0], alpha, L_star, phi_star, nbar))

    return  nbar                                ##  [(Mpc / h)^-3] 

def lya_logmeanlum(z, logLmin=None, printit=False):
    alpha, L_star, phi_star = Sobral_LyaSch(z)

    if logLmin == None:
      ##  Set lower limit on luminosity integral of L*.                                                                                                                                                                                  
      Lmin   =  L_star

    else:
      Lmin   =  10. ** logLmin

    wmin     =  np.log10(Lmin / L_star)

    dw       =  0.1
    ws       =  np.arange(wmin, wmin + 4., dw)   ##  w = log10(L / L*)
    xs       =  np.log(10) * ws                  ##  x =    ln(L / L*)
    ys       =  10. ** ws                        ##  y =       L / L*  

    nbar     =  lya_nbar(z, logLmin=logLmin, printit=printit)  ##  [(Mpc / h)^-3] 

    meanlum  =  phi_star * np.exp(-ys) * 10.**((1. + alpha) * ws) * ys * L_star
    meanlum  =  np.log(10.) * np.sum(meanlum) * dw
    meanlum /=  nbar
    meanlum  =  np.log10(meanlum)
    
    if printit:
      print('%.3lf  %.3lf  %.3lf  %.3le  %.3le  %.3le' % (z, synfits[ind,0], alpha, L_star, phi_star, meanlum))

    return  meanlum  ##  [ergs / s]   
   

if __name__ == '__main__':
    print('\n\nWelcome to a Lyman-alpha Schecter calculator.\n\n')

    logLya  =  42.5

    zs      =  np.arange(2.2, 5.4,  0.1)
    dVs     =  dVols(zs, cosmo, params)

    results =  [] 
    
    for z in zs:
      nbar       = lya_nbar(z, logLmin=logLya, printit=False)      ## [(Mpc / h)^-3]   
      logmeanlum = lya_logmeanlum(z, logLmin=None, printit=False)  ## [ergs / s]
      logmhalo   = logMhalo(logLya, z)                             ## [Msun / h]

      results.append([np.log10(nbar), logmeanlum, logmhalo])

    results = np.array(results)
      
    pl.plot(zs, results[:,0], label=r'$\bar n(z) \ [(h^{-1} \rm{Mpc})^3]$ for $\rm{log}_{10}|L_{\rm{min}}| =$' + ' %.2lf' % logLya)
    ##  pl.plot(zs, results[:,1], label='log$_{10}$|<L> / (ergs/s)|')
    ##  pl.plot(zs, results[:,2], label=r'log$_{10}|M_{\rm{halo}} / (M_\odot / h$)|') 

    pl.xlabel(r'$z$')
    ## pl.ylabel(r'$\bar n(z, L > L_*)$')
    pl.yscale('linear')
    pl.legend()
    pl.show()
    
    print('\n\nDone.\n\n')
