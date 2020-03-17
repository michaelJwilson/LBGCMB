import  numpy                    as      np
import  astropy.units            as      u

from    cosmo                    import  cosmo
from    params                   import  get_params
from    schechterfn              import  SchechterLfn, SchechterMfn
from    dVols                    import  dVols


params = get_params()

def visibilecut(x, xlim, type='mag'):
  '''                                                                                                                                                      
  Function to return unity if M <= M_lim or L >= L_lim; Otherwise, zero.                                                                                   
  '''

  if type == 'mag':
    return  np.piecewise(x, [x >  xlim, x <= xlim], [0., 1.])

  elif type == 'lum':
    return  np.piecewise(x, [x >= xlim, x <  xlim], [1., 0.])

  else:
    raise  ValueError("Requested type is not available")

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
                                                                   ## Gives an M_AB at 1700 \AA for the mean redshift z \tilde 3.05                     
                                                                   ## if NOT k-corrected; otherwise M_AB(v_obs) i.e. median frequency of R.             

  if M_standard == None:
    return  Mlim

  else:
    Llim  = 10.0 ** (-0.4 * (Mlim - M_standard))                   ## Units are luminosity equivalent of 'M_standard'; e.g. M_* gives [L_*].        
                                                                   ## Lv dv = 4 pi D_L^2 F_obs dv_obs 10**(-m/2.5); similary for L* and m*.              
                                                                   ## Gives (Lv/L) = 10**-0.4(m - m*) = 10**-0.4(M - M*) 
    ##  Rest absolute magnitude.
    return  Mlim, Llim

def comovdensity(z, phi_star, M_star, alpha, type='app', mlim=25.0, band='g', printit=True):
  ''' 
  Given a redshift, Schecter fn. parameters and a selection function type, 
  e.g. apparent mag. limited or Goldrush selection, return the expected 
  number density of observed objects.  
  '''

  if type == 'qso':
    from  luminosity_fn  import  gmag
    from  luminosity_fn  import  get_ns


    dM     = 0.05
    Ms     = np.arange(-32.,   -15., dM)
    gs     = gmag(Ms + dM, z, restM=False, printit=False)

    nbar   = get_ns(Ms, zee=z)
    nbar   = np.log10(nbar)

    ##  Cut to maglim.                                                                                                                              
    nbar   = nbar[gs <= mlim][-1]

    return  nbar
  
  elif type == 'lya':
    from  sobral  import  lya_nbar

    ##  Treat mlim as logLmin limit. i.e. Set lower limit 
    ##  on luminosity integral of L* if logLmin == None. 
    ##  Else logLmin.  Return:  ##  [(Mpc / h)^-3]

    result = lya_nbar(z, logLmin=mlim, printit=False)
    result = np.log10(result)

    if printit:
      print(mlim, result)

    return  result

  elif type == 'oii':
    ##  Treat mlim as logLmin [ergs/s] limit.
    from  comparat  import  oii_nbar

    result = oii_nbar(z, logLmin=mlim, printit=False)
    result = np.log10(result)

    if printit:
      print(mlim, result)

    return  result

  else:
    ##  Derived from a Schechter fn.
    import  astropy.constants  as      const
    from    schechterfn        import  SchechterMfn


    ##  Rest-frame absolute magnitude, z=0 and D_L = 10pc.
    MM                 = np.linspace(M_star - 15., M_star + 15., 1000)     ## Integrate from sources 10**10. times brighter than M* to 10**10. times dimmer.
                                                                           ## Even spacing in Magnitude -> intergral dM.     

    PhiMUV             = SchechterMfn(MM, phi_star, M_star, alpha)         ## Note:  Schechter fn. integrals are the incomplete gamma fn.
    
    Mlim, Llim         = mlimitedM(z, mlim, M_star, kcorr=True)            ## Returns [L_standard], e.g. M_* gives [L_*].
    
    if type   == 'app':
      PhiMUV          *= visibilecut(MM, Mlim, type='mag')

    '''
    elif type == 'hsc': 
      hsc_selectiondict = colourcuts.get_SubaruSelectionfn(printit = False, plotit = False) 
      hsc_selection     = hsc_selectiondict[band]
      
      ## Convert each MM (rest-frame absolute AB magnitude) to L_uv(\nu) [erg/s/Hz] @ 2300 A.
      MAB_zero          = -2.5 * np.log10(3631.e-23)

      Lv                = 4. * np.pi * (10. * const.pc.value * 100.)**2. * 10.**(-0.4 * (MM + MAB_zero))
      Lv                = np.log10(Lv)                                      ## np.log10(L_uv(\nu)) [erg/s/Hz] @ 1500 A (2300 A).

      for i, x in enumerate(PhiMUV): 
        PhiMUV[i]      *= np.maximum(hsc_selection((z, Lv[i])), 0.0)        ## interpolation does not enforce positive definite selection.
    '''

    ##  Calculate expected number density (integral dM). 
    nbar     = np.trapz(PhiMUV, MM)                    ## log_10(nbar [(h_100/Mpc)^3]) 
    nbar     = np.log10(nbar)
    
    if printit:
      print("mlim:  %3.3f \t Mlim:  %3.3f \t Llim:  %3.3f \t log10|<n>|:  %3.3f" % (mlim, Mlim, Llim, nbar))
    
    return  nbar

def dndz(zs, phi_star, M_star, alpha, mlim, type='app', printit = True, completeness=None, app_linelim=False):  
  zs, dVs  =  dVols(zs, cosmo, params)
  ns       =  []

  for i, zee in enumerate(zs):
    '''
    Get expected number density for app. mag. (dropout colour) selected expected sample at this redshift slice. 
    Note:  Neglects evolution in luminosity fn., but includes change in app. mag with z.
    '''
    if app_linelim:
      ##  mlim interpreted as a line flux density limit in ergs/s/cm2.
      ##  Note:  lack of h!  Matches DESI science report (on ELGs.)
      ##  eqn. (3.89) of Cosmological Phyiscs: line-flux emission.                                                                                          
      lumdist  = (1. + zee) * cosmo.comoving_distance(zee).value       ##  [Mpc]                                                                        
      lumdist *= 1.e6 * u.parsec.to('cm')                              ##   [cm]                                                                         

      limit    = mlim + np.log10(4. * np.pi) + 2. * np.log10(lumdist)

      nbar     = comovdensity(zee, phi_star, M_star, alpha, type=type, mlim=limit, printit=False)

    else:
      nbar     = comovdensity(zee, phi_star, M_star, alpha, type=type, mlim=mlim, printit=False)  

    nbar = 10. ** nbar               ##  [(h_100/Mpc)^3]

    ns.append(nbar)

  ns = np.array(ns)

  return  zs, dVs, ns

def projdensity(zmin, zmax, phi_star, M_star, alpha, mlim, type='app', printit = True, completeness=None, app_linelim=False):
  '''                                                                                                                                                     
  Integrate (e.g. app. mag. selected) \bar n(z) over a redshift slice of                                                                                  
  width dz to get expected galaxies per sq. degree.                                                                                                       
  '''

  pnbar       = 0.0

  zs          = np.linspace(zmin, zmax, 1500)
  zs, dVs, ns = dndz(zs, phi_star, M_star, alpha, mlim, type=type, printit=printit, completeness=completeness, app_linelim=app_linelim)

  for ii, zee in enumerate(zs):
    if completeness is not None:
      pnbar  += dVs[ii] * ns[ii] * completeness(zee)

    else:
      pnbar  += dVs[ii] * ns[ii]   ##  Aggregate number of each slice.                                                                                  

  pnbar /= 4. * np.pi              ##  Galaxies per steradian.                                                                                          
  pnbar /= (180. / np.pi) ** 2.    ##  Galaxies per sq. degree.                                                                                            

  if printit:
    print('mlim:  %3.3lf \t z:  %.1lf \t %6.6le g/deg^2 \t\t Vol.:  %6.3le (Mpc/h)^3' % (mlim, (zmin + zmax) / 2., pnbar, dVs.sum()))

  return  pnbar
  

if __name__ == "__main__":
  import  pylab         as      pl

  from    reddy.specs   import  samplestats as rsamplestats
  from    specs         import  samplestats as gsample_stats
  from    Malkan.specs  import  samplestats as usample_stats


  print('\n\nWelcome to a Schechter fn. calculator for the projected density of LBG dropouts.\n\n')
                                                                                                 
  dz               =    0.9  
  mlim             =  25.00                          

  for name, key, stats in zip(['BX', 'u', 'g'], ['BX', 'Malkan', 'g'], [rsamplestats(), usample_stats(), gsample_stats()]):    
    print('\n\n------------------------  {}  ------------------------'.format(key))

    midz           =  stats[key]['z']
    alpha          =  stats[key]['schechter']['alpha']
    M_star         =  stats[key]['schechter']['M_star']
    phi_star       =  stats[key]['schechter']['phi_star']

    mlims          =  np.arange(22.5, 26.0, 0.05)
    nbars          =  []
    
    for mlim in mlims:
      nbar         =  comovdensity(midz, phi_star, M_star, alpha, type='app', mlim=mlim, printit=True)
      nbars.append(nbar)

    nbars          =  np.array(nbars)

    pl.plot(mlims, nbars, label=key)

    np.savetxt('simba/dat/nbar_{}.txt'.format(key), np.c_[mlims, nbars], fmt='%.6lf')
    
  pl.legend(frameon=False)

  pl.savefig('simba/plots/simba_nbar.pdf')
  
  print('\n\nDone.\n\n')
