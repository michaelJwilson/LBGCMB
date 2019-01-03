import  os
import  numpy              as      np
import  matplotlib         as      mpl
import  matplotlib.pyplot  as      plt

from    utils              import  comoving_distance
from    pmh                import  Pmm
from    scipy.integrate    import  simps
from    lensing            import  lensing_kernel, Ckk
from    params             import  get_params
from    Cgg                import  Cgg, Ngg
from    Nkk                import  Nkk
from    sliced_pz          import  sliced_pz
from    numpy.linalg       import  inv
from    matplotlib.mlab    import  bivariate_normal
from    fisher_contour     import  plot_ellipse
from    utils              import  latexify 
from    get_dndz           import  get_uvudf_dndz


if __name__ == "__main__":
  import  pylab              as      pl

  from    prep_Llls          import  prep_Llls
  from    Gaussian_pz        import  Gaussian_pz
  from    snr                import  snr
  from    pmh                import  get_PkInterps, linz_bz
  from    prep_camb          import  CAMB, Clxy
  from    bolometers         import  bolometers
  from    completeness       import  get_dropoutpz

  from    Malkan.specs       import  samplestats         as  usample_stats
  from    specs              import  samplestats         as  gsample_stats
  from    ilim               import  get_nbar_nocontam
  from    schmittfull        import  ss_pz
  from    scipy.interpolate  import  interp1d
  from    bz                 import  get_dropoutbz
  from    Ckg                import  Ckg, var_Ckg
  from    planck18_bao       import  get_sig8z


  print("\n\nWelcome to a calculator for the (parameter) bias due to the interloper dN/dz.\n\n")
  
  ##  Set CMB configuration. 
  cmbexp                             = 'CMBS4'

  fsky, thetab, DeltaT, iterative    =  bolometers[cmbexp]['fsky'],   bolometers[cmbexp]['thetab'],\
                                        bolometers[cmbexp]['DeltaT'], bolometers[cmbexp]['iterative']
  
  cambx                              =  CAMB()

  NLlls, Llls, nmodes                =  prep_Llls(NLlls = 60, Lmin = 50., Lmax = 5000., log10=True)
  Pk_interps                         =  get_PkInterps(cambx)
  
  ##  No Detector noise -- this is handled by Clxy of prep_camb.                                                                                    
  (lensCl_interps, nolensCl_interps) = cambx.get_Cls()
  
  ##  Galaxies per sq. degree.                                                                                                                        
  band                 =  'u'

  if band == 'g':
    stats              =  gsample_stats()
    stats              =  get_nbar_nocontam(band, depth='W', printit=False)

    nbar               =  stats[band]['nbar_nointerlopers']
    peakz              =  stats[band]['z']

  elif band == 'u':
    stats              =  usample_stats()
    ## stats           =  get_nbar_nocontam(band, depth='W', printit=False)

    ## nbar            =  stats['Malkan']['nbar']
    peakz              =  stats['Malkan']['z']
    
    ##  Get dropout Schechter counts for given band.                                                                                                         
    root               =  os.environ['LBGCMB']
    data               =  np.loadtxt(root + "/dropouts/nz/schechter/dat/schechter_estimate_Malkan_dropouts.txt")

    ##  Angular density on the sky with mag. depth. 
    ms                 =  data[:,0][::-1]
    Npz                =  data[:,1][::-1]

    nbar               =  Npz[0]

  else:
    raise ValueError('\n\nChosen band is not available.\n\n')

  zmin                 =  peakz - 2.00
  zmax                 =  peakz + 2.00

  drop_bz              =  get_dropoutbz()
  bz                   =  lambda z:  drop_bz(peakz)

  ##  Get Rafelski based estimates on p(z).
  pz                   =  get_uvudf_dndz(band='g', field='UVUDF', mag='')

  ##  Generate Ckg. 
  ckg                  =  Ckg(Pk_interps, Llls, zmin, zmax, pz, bz, zeff=False)
  nkk                  =  Nkk(lensCl_interps, nolensCl_interps, Llls, terms=['TT', 'TE', 'EE', 'EB'],\
                              thetab=thetab, DeltaT=DeltaT, iterative=iterative, pickleit=False)

  vkg                  =  var_Ckg(Pk_interps, lensCl_interps, nolensCl_interps, Llls, zmin, zmax, pz, bz, nbar, fsky,\
                                  nkk=nkk, iterative=iterative, pickleit=True)

  snr                  =  snr(Llls, ckg, vkg, nmodes, lmax=1.e5)

  print('\n\nTotal S/N for Ckg:  %.1lf\n\n' % snr)
  
  zs, sig8z, esig8z    = get_sig8z(interp=True)

  ##  and the parameter Fisher matrix.
  fid_b1               =  drop_bz(peakz)
  fid_sig8             =    sig8z(peakz)
  
  ##  Inversion check                                                                                                                                     
  ##  Fisher           =  np.array([[4., 7.], [2., 6.]])                                                                                                  
  
  Fisher               =  np.array([[4. / fid_sig8 ** 2., 2. / fid_b1 / fid_sig8], [2. / fid_b1 / fid_sig8, 1. / fid_b1 ** 2.]])
  Fisher              *=  snr ** 2.

  iFisher              =  inv(Fisher)

  print("\n\nFisher:\n")
  print   Fisher
  
  print("\n\nand its inverse:\n")
  print  iFisher

  ##  Change in Ckg with dN/dz -> dN/dz' 
  bzz        =  bz
  pzz        =  get_uvudf_dndz(band='g', field='UVUDF', mag='24.5')

  dckg       =  ckg - Ckg(Pk_interps, Llls, zmin, zmax, pzz, bzz, zeff=True)

  print("\n\n(dCkg / Ckg): \n\n%s \n\nper cent" % np.array_str(100. * dckg / ckg))

  ##  C^{-1} dCkg
  ndckg      = dckg / vkg 
  
  vec        = np.array([2. / fid_sig8, 1. / fid_b1]) * np.sum(ckg * ndckg)  
  dtheta     = np.dot(iFisher, vec)

  print("\n\nd(sig8): %.6lf  and  d(b1): %.6lf" % (dtheta[0], dtheta[1]))

  
  ##  And plot ... 
  latexify(fig_width=None, fig_height=None, columns=1, equal=True)

  fig = plt.gcf()
  ax  = fig.add_subplot(111)
  
  pl.plot(fid_b1, fid_sig8, 'k*', label='Fiducial', markersize=5)

  plot_ellipse(x_cent = fid_b1 + dtheta[1], y_cent = fid_sig8 + dtheta[0], ax = ax, cov = iFisher, mass_level = 0.99,\
               fill=False, fill_kwargs={'alpha':  0.5}, plot_kwargs={'c': 'g'})

  plot_ellipse(x_cent = fid_b1 + dtheta[1], y_cent = fid_sig8 + dtheta[0], ax = ax, cov = iFisher, mass_level = 0.95,\
               fill=False, fill_kwargs={'alpha':  0.5}, plot_kwargs={'c': 'r'})

  plot_ellipse(x_cent = fid_b1 + dtheta[1], y_cent = fid_sig8 + dtheta[0], ax = ax, cov = iFisher, mass_level = 0.68,\
               fill=False, fill_kwargs={'alpha':  0.5}, plot_kwargs={'c': 'b'})
    
  ## plt.arrow(fid_b1, fid_sig8, dtheta[1], dtheta[0], 'k', linewidth=.25, head_width=0.01, overhang=.8, color='k')

  pl.xlim(0.95 * fid_b1,   1.05 * fid_b1)
  pl.ylim(0.95 * fid_sig8, 1.05 * fid_sig8)

  pl.xlabel(r'$b_1$')
  pl.ylabel(r'$\sigma_8$')

  pl.savefig('plots/interloper_bias_%sband.pdf' % band, bbox_inches='tight')
  
  print("\n\nDone.\n\n")
