import  matplotlib; matplotlib.use('PDF')

import  numpy                    as      np
import  astropy.units            as      u
import  pylab                    as      pl
import  pandas                   as      pd

from    astropy.table            import  Table
from    reddy.specs              import  samplestats as rsamplestats
from    specs                    import  samplestats as gsample_stats
from    Malkan.specs             import  samplestats as usample_stats
from    cosmo                    import  cosmo
from    params                   import  get_params
from    schechterfn              import  SchechterLfn, SchechterMfn
from    dVols                    import  dVols


params = get_params()

def  mlimitedM(z, mlim, M_standard=None, kcorr=True):
  '''                                                                                                                                                                                                                                 
  Return M_lim (L_lim) in units of M_standard (L_standard) for given redshift                                                                                                                                                         
  and apparent mag. limit.  Here, M_standard is e.g. M* for the Schechter fn.                                                                                                                                                         
  '''

  from  utils  import  comoving_distance


  z     = np.asarray(z)
  aa    = 1. / (1. + z)                                            ## scale factor at z.                                                                                                                                               
  chi   = comoving_distance(z)                                     ## [Mpc/h]                                                                                                                                                          
  dmod  = 25. + 5. * np.log10(chi / aa / params['h_100'])          ## \mu = 5.*log_10(D_L/10pc) = 25. + 5 log_10(D_L) for [D_L] = Mpc.                                                                                                
                                                                   ## D_L = (1. + z)*chi.                                                                                                                                             
  if kcorr:
    kcorr =  -2.5 * np.log10(1.0 / aa)                             ## Eq. (14) of Reddy++ and preceding text; assumes flat Fv.
    
  else:
    kcorr =   0.0                                                  ## Assumed zero k-correction.                                                                                                                                       
  Mlim  = mlim - dmod - kcorr                                      ## Apparent mag. limited (~5 sig. source detected) in R.
                                                                   ## Gives an M_AB at 1700 \AA for the mean redshift z \tilde 3.05                                                                                                                                                                      ## if NOT k-corrected; otherwise M_AB(v_obs) i.e. median frequency of R.                                                                                           
  if M_standard == None:
    return  Mlim

  else:
    Llim  = 10.0 ** (-0.4 * (Mlim - M_standard))                   ## Units are luminosity equivalent of 'M_standard'; e.g. M_* gives [L_*].                                                                                          
                                                                   ## Lv dv = 4 pi D_L^2 F_obs dv_obs 10**(-m/2.5); similary for L* and m*.                                                                                           
                                                                   ## Gives (Lv/L) = 10**-0.4(m - m*) = 10**-0.4(M - M*)                                                                                                              
    ##  Rest absolute magnitude.                                                                                                                                                                                                      
    return  Mlim, Llim

def PhiMUV(z, phi_star, M_star, alpha):
  ''' 
  Given a redshift, Schecter fn. parameters and a selection function type, 
  e.g. apparent mag. limited or Goldrush selection, return the expected 
  number density of observed objects.  
  '''

  ##  Derived from a Schechter fn.
  import  astropy.constants  as      const
  from    schechterfn        import  SchechterMfn
    

  ##  Rest-frame absolute magnitude, z=0 and D_L = 10pc.
  MM                 = np.linspace(M_star - 15., M_star + 15., 1000)     ## Integrate from sources 10**10. times brighter than M* to 10**10. times dimmer.

  ##  Even spacing in Magnitude -> intergral dM.     
  PhiMUV             = SchechterMfn(MM, phi_star, M_star, alpha)         ## Note:  Schechter fn. integrals are the incomplete gamma fn.
 
  Mlims              = np.arange(-22., -16., 0.1)
  result             = np.array([np.trapz(PhiMUV[MM < Mlim], MM[MM < Mlim]) for Mlim in Mlims])
  
  return  Mlims, result


if __name__ == "__main__":
  print('\n\nWelcome to a Schechter fn. calculator for the projected density of LBG dropouts.\n\n')
                                                                                                 
  dz               =    0.9  

  nodust           =  pd.read_csv('simba/dat/no_dust_Ns.csv')
  dust             =  pd.read_csv('simba/dat/dust_Ns.csv')

  output           =  Table()
  
  for z, mlim, c, name, key, stats in zip([2.0, 3.0, 4.0], [25.5, 24.6, 25.8, 25.8], ['gold', 'b', 'g'], ['BX', 'u', 'g'], ['BX', 'Malkan', 'g'], [rsamplestats(), usample_stats(), gsample_stats()]):    
    print('\n\n------------------------  {}  ------------------------'.format(key))

    midz           =  stats[key]['z']
    alpha          =  stats[key]['schechter']['alpha']
    M_star         =  stats[key]['schechter']['M_star']
    phi_star       =  stats[key]['schechter']['phi_star']
    
    Mlims, result  =  PhiMUV(midz, phi_star, M_star, alpha)

    pl.semilogy(Mlims, result, label=r'${}$'.format(key), c=c, alpha=0.3)
    
    output[key + '_MLIM']  = Mlims
    output[key]            = result
    
    print(key, phi_star, M_star, alpha, mlimitedM(z, mlim, M_standard=None, kcorr=True))

    output.meta[key[:4].upper() + '_ML'] = mlimitedM(z, mlim, M_standard=None, kcorr=True)
    
  output.pprint(max_width=-1)

  output.write('cumlumuv.fits', format='fits', overwrite=True)
  
  pl.legend(frameon=False, ncol=3, loc=4)
  
  pl.savefig('simba/plots/simba_cumlumuv.pdf')
  
  print('\n\nDone.\n\n')

  
