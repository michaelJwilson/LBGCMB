import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    Cgg                import  Cgg, Ngg
from    Ckg                import  Ckg, var_Ckg
from    prep_Llls          import  prep_Llls
from    Gaussian_pz        import  Gaussian_pz
from    pmh                import  Pmm, get_PkInterps, linz_bz
from    prep_camb          import  CAMB
from    schmittfull_nz     import  ss_pz, load_samples
from    schmittfull_bz     import  get_allbz
from    completeness       import  get_dropoutpz
from    specs              import  samplestats                                                                                                             
from    ilim               import  get_nbar_nocontam
from    scipy.interpolate  import  interp1d
from    prep_camb          import  CAMB
from    lensing            import  Ckk, var_Ckk
from    Nkk                import  Nkk
from    bolometers         import  bolometers
from    zeldovich_Lmax     import  Lcutmax
from    pylab              import  rcParams


if __name__ == '__main__':
    print('\n\nWelcome to the cross-correlation calculator.\n\n')
    
    ##  Prepare pycamb module; linear, non-linear matter P(k) and Cls.                                                                                      
    cambx                           =  CAMB()
    Pk_interps                      =  get_PkInterps(cambx)

    NLlls, Llls, nmodes             =  prep_Llls(NLlls = 60, Lmin = 50., Lmax = 5000., log10=True)

    ##  No Detector noise -- this is handled by Clxy of prep_camb.                                                                                          
    (lensCl_interps, nolensCl_interps) = cambx.get_Cls()

    cmbexp                          = 'CMBS4'
    fsky, thetab, DeltaT, iterative =  bolometers[cmbexp]['fsky'], bolometers[cmbexp]['thetab'],\
                                       bolometers[cmbexp]['DeltaT'], bolometers[cmbexp]['iterative']

    ckk                             =  Ckk(Pk_interps, Llls)

    nkk                             =  Nkk(lensCl_interps, nolensCl_interps, Llls, terms=['TT', 'TE', 'EE', 'EB'],\
                                           thetab=thetab, DeltaT=DeltaT, iterative=iterative, pickleit=False)
    
    ## --  Galaxy samples  -- ##                                                                                        
    '''
    band                 =  'g'                                                                                                                           
    stats                =  samplestats()                                                                                                           
    stats                =  get_nbar_nocontam(band, depth='W', printit=False)                                                                              

    peakz                =  stats[band]['z']                                                                                                         
    nbar                 =  stats[band]['nbar_nointerlopers']                                                                                              
    
    
    ##  Effectively overwrites hard z limits above.                                                                                                         
    zee, pzee            =  get_dropoutpz()                                                                                                               
    pz                   =  interp1d(zee, pzee, kind='linear', bounds_error=False, fill_value=0.0, assume_sorted=False)
    
    
    zmin                 =  peakz - 2.00
    zmax                 =  peakz + 2.00
    '''
    
    ## --  Schmittfull and Seljak  -- ##                                                                                
    df = load_samples()
    
    ##  Ordered dictionary with keys corresponding to row['Surveys'].
    bs = get_allbz()
    
    for index, row in df.iterrows():
      zmin    =  row['zmin']
      zmax    =  row['zmax']     
      nbar    =  row['N']

      print('\t%s \t %.2lf \t %.2lf \t %.6le' % (row['Surveys'], zmin, zmax, nbar))
      
      ##  Assumed uniform between zmin and zmax. 
      peakz   =  0.5 * (zmin + zmax)
        
      dz      =  0.01
      zs      =  np.arange(0.0, 5.0, dz)

      result  =  np.zeros_like(zs)
      result[(zs >= zmin) & (zs <= zmax)] = 1.0
        
      result /=  (np.sum(result) * dz) 
      pz      =  interp1d(zs, result, kind='linear', bounds_error=False, fill_value=0.0, assume_sorted=False)

      ckg     =  Ckg(Pk_interps, Llls, zmin, zmax, pz, bs[row['Surveys']], zeff=True)
      cgg     =  Cgg(Pk_interps, Llls, zmin, zmax, pz, bs[row['Surveys']], zeff=True)

      ngg     =  Ngg(Llls, zmin, zmax, pz, nbar)
      
      rho     = ckg / np.sqrt((cgg + ngg) * (ckk + nkk))

      pl.plot(Llls, rho, label = row['Surveys'])

    pl.xlabel(r'$L$')
    pl.ylabel(r'$\rho_L$')

    pl.ylim(0.0, 1.0)

    pl.legend(ncol=1, loc=1)

    plt.tight_layout()

    pl.savefig('plots/rho.pdf')
    
    print('\n\nDone.\n\n')
