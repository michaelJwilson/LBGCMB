import  os
import  matplotlib         as      mpl
import  matplotlib.pyplot  as      plt
import  numpy              as      np
import  astropy.units      as      u

from    cosmo              import  cosmo
from    params             import  get_params
from    utils              import  comoving_distance
from    astropy            import  constants as const
from    scipy.integrate    import  simps
from    pmh                import  Pmm
from    dropouts.reddy     import  samplestats
from    sliced_pz          import  sliced_pz
from    utils              import  latexify


latexify(fig_width=None, fig_height=None, columns=1, equal=True, fontsize=10)

params   = get_params()

def snr(Llls, Cls, Vars, nmodes, lmax=1.e5):
  result  = (Cls[Llls <= lmax] / np.sqrt(Vars[Llls <= lmax])) ** 2.
  result  = np.sum(nmodes[Llls <= lmax] * result)

  return  np.sqrt(result)


if __name__ == '__main__':
  import  pylab              as      pl

  from    prep_Llls          import  prep_Llls
  from    Gaussian_pz        import  Gaussian_pz
  from    pmh                import  Pmm, get_PkInterps, linz_bz
  from    prep_camb          import  CAMB
  ##  from    schmittfull        import  ss_pz
  from    completeness       import  get_dropoutpz
  
  ## from    specs              import  samplestats
  from    ilim               import  get_nbar_nocontam
  from    scipy.interpolate  import  interp1d
  from    prep_camb          import  CAMB
  from    bolometers         import  bolometers
  from    lensing            import  Ckk, var_Ckk
  from    Nkk                import  Nkk
  from    Ckg                import  Ckg, var_Ckg
  from    bolometers         import  bolometers
  from    zeldovich_Lmax     import  Lcutmax
  from    pylab              import  rcParams
  from    bz                 import  get_dropoutbz


  print('\n\nWelcome to snr.\n\n')
  '''
  ##  Galaxies per sq. degree.                                                                                                                            
  band                 =  'g'
  stats                =  samplestats()
  stats                =  get_nbar_nocontam(band, depth='W', printit=False)

  peakz                =  stats[band]['z']
  nbar                 =  stats[band]['nbar_noint']                                                                                               

  ##  Effectively overwrites hard z limits above.                                                                                                          
  zee, pzee            =  get_dropoutpz()                                                                                                                   
  pz                   =  interp1d(zee, pzee, kind='linear', bounds_error=False, fill_value=0.0, assume_sorted=False)                                   
  '''
  
  ##  Reddy u-drops.
  band                 =  'LBG'
  stats                =  samplestats()
  
  peakz                =  stats[band]['z']
  nbar                 =  stats[band]['nbar_noint'] 
  
  ##  Defaults to Hildebrandt (2009).
  pz                   =  Gaussian_pz
  
  ##  Stored z, b(z) for LBGs at z=3, 4 etc. (u and g respectively).
  bz                   =  linz_bz 
  ##  bz               =  get_dropoutbz()
  
  ##  Schmittfull and Seljak (2017).                                                                                                                          
  ##  pz, nbar         =  ss_pz()

  zmin                 =  peakz - 2.00
  zmax                 =  peakz + 2.00

  Lmax                 =  Lcutmax[np.round(peakz)][0] 

  ##  Prepare pycamb module; linear, non-linear matter P(k) and Cls.                                                                        
  cambx                =  CAMB()
  NLlls, Llls, nmodes  =  prep_Llls(NLlls = 60, Lmin = 50., Lmax = 5000., log10=True)

  Pk_interps           =  get_PkInterps(cambx)

  ## No Detector noise -- this is handled by Clxy of prep_camb.                                                                                   
  (lensCl_interps, nolensCl_interps) = cambx.get_Cls()

  ## Dilution factors.
  dfactors                           =  np.array([0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 1.e1, 5.e1])
  
  for cmbexp in bolometers: 
    results = []

    ckg     = Ckg(Pk_interps, Llls, zmin, zmax, pz, bz, zeff=True)

    fsky, thetab, DeltaT, iterative = bolometers[cmbexp]['fsky'],   bolometers[cmbexp]['thetab'],\
                                      bolometers[cmbexp]['DeltaT'], bolometers[cmbexp]['iterative']

    nkk     = Nkk(lensCl_interps, nolensCl_interps, Llls, terms=['TT', 'TE', 'EE', 'EB'],\
                  thetab=thetab,  DeltaT=DeltaT,    iterative=iterative, pickleit=False)

    for dfactor in dfactors:    
      dnbar               = dfactor * nbar
  
      vkg                 =  var_Ckg(Pk_interps, lensCl_interps, nolensCl_interps, Llls, zmin, zmax, pz, bz, dnbar, fsky,\
                                     nkk=nkk,    iterative=iterative, pickleit=True)
      
      result              =  snr(Llls, ckg, vkg, nmodes, lmax = Lmax)
      result             /=  np.sqrt(fsky)
      
      results.append(result)

      print('\n\nFor z=%.1lf, total S/N: %.1lf to Lmax of %.1lf (fsky, thetab, DeltaT = %.3lf, %.2lf, %.2lf)' % (peakz, result, Lmax, fsky, thetab, DeltaT))

    results = np.array(results) 
    
    pl.semilogx(dfactors * nbar, results, label=cmbexp)
  
  pl.xlabel(r'$\bar n / \rm{deg}^{2}$', fontsize=12)
  pl.ylabel('Cumulative (S/N) ' + r'$/ \ \sqrt{f_{\rm{sky}}}$', fontsize=12)

  pl.xlim(1.e1, 1.e4)
  pl.ylim(0.0,  550.)
  
  pl.legend(loc=2, ncol=2)

  plt.tight_layout()
  
  pl.savefig('plots/%ssnr.pdf' % band)

  print('\n\nDone.\n\n')

