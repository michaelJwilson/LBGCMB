import  os
import  matplotlib          as      mpl
import  matplotlib.pyplot   as      plt
import  numpy               as      np
import  astropy.units       as      u
  
from    params              import  get_params
from    pmh                 import  get_PkInterps, linz_bz
from    bz                  import  get_dropoutbz
from    utils               import  latexify
from    nbar                import  projdensity
from    zeldovich_Lmax      import  Lcutmax
from    Cgg                 import  Cgg, var_Cgg, Ngg
from    numpy.linalg        import  inv
from    fisher_contour      import  plot_ellipse
from    get_uvudf_pz        import  get_uvudf_pz
from    planck18_bao        import  get_sig8z
from    whitebook_pz        import  whitebook_pz


params = get_params()

def get_allCls(Pk_interps, Llls, nbar, fsky, zmin, zmax, pz, bz, bz2 = None, survey_pz2 = None, zeff=True, samplevar_lim=False):
  ##  Nkk.
  nkk          =      Nkk(lensCl_interps, nolensCl_interps, Llls, terms=['TT', 'TE', 'EE', 'EB'],\
                          thetab=thetab, DeltaT=DeltaT, iterative=iterative, pickleit=False)
  ##  Cgg.                                    
  cgg          =      Cgg(Pk_interps, Llls, zmin, zmax, pz, bz, zeff=zeff, bz2 = bz, survey_pz2 = pz)
  ngg          =      Ngg(Llls, zmin, zmax, pz, nbar)
  vgg          =  var_Cgg(Llls, zmin, zmax, pz, bz, nbar, fsky, samplevar_lim=samplevar_lim, cgg = cgg)

  ##  Ckg.                                                                                                                                                  
  ckg          =      Ckg(Pk_interps, Llls, zmin, zmax, pz, bz, zeff=zeff)
  vkg          =  var_Ckg(Pk_interps, lensCl_interps, nolensCl_interps, Llls, zmin, zmax, pz, bz, nbar, fsky,\
                          nkk=nkk, iterative=iterative, pickleit=True)
  ##  Ckk                                                                                                                                                   
  ckk          =      Ckk(Pk_interps, Llls)
  vkk          =  var_Ckk(Llls, fsky, nkk, Pk_interps, samplevar_lim=samplevar_lim)
  
  return  {'gg': cgg + ngg, 'kg': ckg, 'gk': ckg, 'kk': ckk + nkk}
  
def get_ClsCov(Llls, Cls, fsky):
  result = {}

  for key in Cls:
    for kkey in Cls:
      result[key + kkey] = (Cls[key[0] + kkey[0]] * Cls[key[1] + kkey[1]] + Cls[key[0] + kkey[1]] * Cls[key[1] + kkey[0]]) / (2. * Llls + 1.) / fsky
  
  for key in result:
    print key, result[key]

  return  result

def dp_DL(p, fid_sig8, fid_b1, DL):
  ##  DL = (Ckk, Ckg, Cgg) for given L.  Partial derivative wrt the parameters {sig8, b1}.
  if   p == 's8':
    return  (2. / fid_sig8) * DL

  elif p == 'b1':
    return   np.array([0.0, DL[1] / fid_b1, 2. * DL[2] / fid_b1])
    
  else:
    raise  ValueError('\n\nDerivative with respect to given parameter is not available.\n\n')


if __name__ == "__main__":
  import  pylab              as      pl

  from    prep_Llls          import  prep_Llls
  from    Gaussian_pz        import  Gaussian_pz
  from    prep_camb          import  CAMB
  from    schmittfull        import  ss_pz
  from    completeness       import  get_dropoutpz       as  get_gdropoutpz
  from    Malkan.specs       import  samplestats         as  usample_stats                                                                                 
  from    specs              import  samplestats         as  gsample_stats
  from    ilim               import  get_nbar_nocontam
  from    scipy.interpolate  import  interp1d
  from    bolometers         import  bolometers
  from    lensing            import  Ckk, var_Ckk
  from    Nkk                import  Nkk
  from    Ckg                import  Ckg, var_Ckg 


  print("\n\nWelcome to a calculator for the parameter bias due to interlopers.\n\n")
  
  ## Prepare pycamb module; linear, non-linear matter P(k) and Cls.                                                                                     
  cambx                              =  CAMB()
  Pk_interps                         =  get_PkInterps(cambx)

  ## No Detector noise -- this should be handled by Clxy of prep_camb.                                                                                  
  (lensCl_interps, nolensCl_interps) =  cambx.get_Cls()

  NLlls, Llls, nmodes                =  prep_Llls(NLlls = 60, Lmin = 50., Lmax = 5000., log10=True)

  cmbexp                             =  'CMBS4'

  fsky, thetab, DeltaT, iterative    =  bolometers[cmbexp]['fsky'],   bolometers[cmbexp]['thetab'],\
                                        bolometers[cmbexp]['DeltaT'], bolometers[cmbexp]['iterative']

  ##  Get interp1d for sigma8(z) and its Planckian error.                                                                                               
  zs, sig8z, esig8z                  =  get_sig8z(interp=True)

  ##  Dropout selection.   
  band = 'g'
  
  if  band == 'g':
    stats        =  gsample_stats()
    stats        =  get_nbar_nocontam(band, depth='W', printit=False)
  
    ##  Effectively overwrites hard z limits above.                                                                                                 
    zee, pzee    =  get_gdropoutpz()
    pz           =  interp1d(zee, pzee, kind='linear', bounds_error=False, fill_value=0.0, assume_sorted=False)

    nbar         =  stats[band]['nbar_nointerlopers']
    peakz        =  stats[band]['z']

    detband      =   'i'                                                                                          ## Detection band.
    colors       =  ['darkgreen', 'limegreen', 'g']
  
  elif band == 'Malkan': 
    ##  Reddy u-drops.
    stats        =  usample_stats()

    peakz        =  stats[band]['z']

    alpha        =  stats[band]['schechter']['alpha']
    Mstar        =  stats[band]['schechter']['M_star']
    phi_star     =  stats[band]['schechter']['phi_star']
    
    dzee         =  0.61 / 2.
    peakz        =  stats[band]['z']
    nbar         =  projdensity(peakz - dzee / 2., peakz + dzee / 2., phi_star, Mstar, alpha, mlim=24.5, printit = True, completeness=None)
  
    pz           =  Gaussian_pz

    decband      =   'R'
    colors       =  ['darkblue', 'deepskyblue', 'b']

  else:
    raise ValueError('\n\nChosen band is not available.\n\n')

  ##  Bias with z.
  drop_bz            =  get_dropoutbz()
  bz                 =  lambda z:  drop_bz(peakz)

  ##  Get Rafelski based estimates on p(z).                                                                                                               
  pz                 =  lambda z:  whitebook_pz(z, ilim = 25.3)
  ## pz              =  get_uvudf_dndz(band='g', field='UVUDF', mag='')    

  ##  Change in Ckg with dN/dz -> dN/dz'.  TO DO:  Set low-z population to have red galaxy bias.                                                                    
  bzz                =  lambda z:  1.0 * bz(z)

  pzz                =  lambda z:  whitebook_pz(z, ilim = 22.3)
  ##  pzz            =  get_uvudf_dndz(band='g', field='UVUDF', mag='24.5') 

  ##  and the parameter Fisher matrix.                                                                                                                    
  fid_b1             =  drop_bz(peakz)
  fid_sig8           =    sig8z(peakz)

  print('\n\nFiducial b1(z) and fsig8(z):  %.3lf and %.3lf.' % (fid_b1, fid_sig8))

  zmin               =   0.01
  zmax               =  10.00 

  ##  {'gg': cgg + ngg, 'kg': ckg, 'gk': ckg, 'kk': ckk + nkk}
  Cls                =  get_allCls(Pk_interps, Llls, nbar, fsky, zmin, zmax, pz, bz,   zeff=False, samplevar_lim=False)

  ## 
  xCls               =  get_allCls(Pk_interps, Llls, nbar, fsky, zmin, zmax, pzz, bzz, zeff=False, samplevar_lim=False)

  pl.semilogy(Llls,  Cls['gg'],  label='base Cgg')
  pl.semilogy(Llls, xCls['gg'],  label='shifted Cgg')

  pl.semilogy(Llls,  Cls['kg'],  label='base Ckg')
  pl.semilogy(Llls, xCls['kg'],  label='shifted Ckg')

  pl.xlim(50., 4.e3)
  pl.xlabel(r'L')

  pl.legend(ncol=2)

  pl.savefig('plots/interloper_bias_cls_%sband.pdf' % band, bbox_inches='tight')  

  pl.clf()


  '''
  ##  'ggkk' etc ...
  covs               =  get_ClsCov(Llls, Cls, fsky)

  ##  Covariance matrix of {kk, kg, gg} for given L.
  Cov_Lll            =  np.zeros(3 * 3).reshape(3,3)
  iCov_Lll           =  np.zeros_like(Cov_Lll) 

  ##  Parameter Fisher matrix.
  Fisher             =  np.zeros(2 * 2).reshape(2, 2)
  iFisher            =  np.zeros_like(Fisher)

  ##  Finally, the parameter bias.                                                                                                                        
  dtheta             =  np.zeros(2)
  interim            =  np.zeros(2)


  for LL, dummy in enumerate(Llls):
    dckg         =  Cls['kg'][LL] - Ckg(Pk_interps, Llls, zmin, zmax, pzz, bzz, zeff=False)[LL]
    dcgg         =  Cls['gg'][LL] - Cgg(Pk_interps, Llls, zmin, zmax, pzz, bzz, zeff=False)[LL] - Ngg(Llls, zmin, zmax, pzz, nbar)[LL]

    ##  Construct D(LL) = [Ckk, Ckg, Cgg].                                                                                  
    DL           =  np.array([Cls['kk'][LL], Cls['kg'][LL], Cls['gg'][LL]])

    for i, tracer in enumerate(['kk', 'kg', 'gg']):
      for j, ttracer in enumerate(['kk', 'kg', 'gg']):
        ##  Restricted to first L=50 mode.
        Cov_Lll[i, j] = covs[tracer + ttracer][LL]

    ## .. and invert.  
    iCov_Lll     =  inv(Cov_Lll)

    for i, pp in enumerate(['s8', 'b1']):
      for j, ss in enumerate(['s8', 'b1']):
        ##  Here, DL restricted to first L=50 mode. 
        Fisher[i, j] = np.dot(dp_DL(pp, fid_sig8, fid_b1, DL), np.dot(iCov_Lll, dp_DL(ss, fid_sig8, fid_b1, DL)))
    
    ##  No change to kk due to interloper population. 
    DDL           =  np.array([0.0, dckg, dcgg])
    XXX           =    np.dot(iCov_Lll, DDL)
    
    for j, pparam in enumerate(['s8', 'b1']):
      interim[j] +=  np.dot(dp_DL(pparam, fid_sig8, fid_b1, DL), XXX)
    ##  End of loop over LL.   
  
  ##  Invert Fisher. 
  iFisher  =  inv(Fisher)
  
  for i, param in enumerate(['s8', 'b1']):
    for j, pparam in enumerate(['s8', 'b1']):
      dtheta[i] += iFisher[i,j] * interim[j]


  ##  And plot contour ...                                                                                                                                        
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

  ## pl.xlim(0.95 * fid_b1,   1.05 * fid_b1)
  ## pl.ylim(0.95 * fid_sig8, 1.05 * fid_sig8)

  pl.xlabel(r'$b_1(z=%.1lf)$'      % peakz)
  pl.ylabel(r'$\sigma_8(z=%.1lf)$' % peakz)

  pl.savefig('plots/interloper_bias_%sband.pdf' % band, bbox_inches='tight')
  '''
  print("\n\nDone.\n\n")
