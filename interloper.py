import  os
import  matplotlib          as      mpl
import  matplotlib.pyplot   as      plt
import  numpy               as      np
import  astropy.units       as      u
  
from    params              import  get_params
from    pmh                 import  get_PkInterps
from    bz                  import  get_dropoutbz
from    utils               import  latexify
from    nbar                import  projdensity
from    zeldovich_Lmax      import  Lcutmax
from    Cgg                 import  Cgg, var_Cgg, Ngg
from    numpy.linalg        import  inv
from    fisher_contour      import  plot_ellipse
from    Raf15               import  Raf15_pz
from    planck18_bao        import  get_sig8z
from    whitebook_pz        import  whitebook_pz


params = get_params()

def get_allCls(Pk_interps, Llls, nbar, fsky, zmin, zmax, pz, bz, bz2 = None, survey_pz2 = None, zeff=False, samplevar_lim=False):
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
  
  return  {'gg': cgg, 'kg': ckg, 'gk': ckg, 'kk': ckk}, {'gg': ngg, 'kg': np.zeros_like(ckg), 'gk': np.zeros_like(ckg), 'kk': nkk}
  
def get_ClsCov(Llls, Cls, Nls, fsky, printit=False):
  result = {}

  for key in Cls:
    for kkey in Cls:
      ##  Gaussian covariance;  eqn. (14) of https://arxiv.org/pdf/1710.09465.pdf.
      result[key + kkey]  = (Cls[key[0] + kkey[0]] + Nls[key[0] + kkey[0]]) * (Cls[key[1] + kkey[1]] + Nls[key[1] + kkey[1]]) 

      ##  Cross-terms should be identically zero.
      result[key + kkey] += (Cls[key[0] + kkey[1]] + Nls[key[0] + kkey[1]]) * (Cls[key[1] + kkey[0]] + Nls[key[1] + kkey[0]])

      ##  To be viewed as covariance sampled at Llls.  I.e. to be corrected
      result[key + kkey] /= (2. * Llls + 1.)                                                                                  

      ##  for further binning in L.
      result[key + kkey] /=  fsky                                                                                            

  if printit:
    for key in result:
      print(key, result[key])

  return  result

def dp_DL(p, fid_sig8, fid_b1, DL):
  ##  DL = (Ckk, Ckg, Cgg) for given L.  Partial derivative wrt the parameters {sig8, b1}.
  if  p == 's8':
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
  
  ##  Prepare pycamb module; linear, non-linear matter P(k) and Cls.                                                                                     
  plotit                             =  True 
  
  cambx                              =  CAMB()
  Pk_interps                         =  get_PkInterps(cambx)
  
  ##  No Detector noise -- this should be handled by Clxy of prep_camb.                                                                                  
  (lensCl_interps, nolensCl_interps) =  cambx.get_Cls()

  NLlls, Llls, nmodes                =  prep_Llls(NLlls = 60, Lmin = 50., Lmax = 5000., log10=True)

  ##  CMB specification.
  cmbexp                             =  'CMBS4'

  fsky, thetab, DeltaT, iterative    =  bolometers[cmbexp]['fsky'],   bolometers[cmbexp]['thetab'],\
                                        bolometers[cmbexp]['DeltaT'], bolometers[cmbexp]['iterative']
  
  ##  Get interp1d for sigma8(z) and its Planckian error.                                                                                               
  zs, sig8z, esig8z                  =  get_sig8z(interp=True)

  ##  Dropout selection.   
  band  =     'u'

  zmin  =   0.01
  zmax  =  10.00
  
  if  band == 'g':
    stats        =  gsample_stats()
    stats        =  get_nbar_nocontam(band, depth='W', printit=False)
  
    peakz        =  stats[band]['z']

    ##  Contamination corrected estimate.
    nbar         =  stats[band]['nbar_noint']
    nbar_wint    =  stats[band]['nbar']

    ##  With and without interlopers;  Sq. deg. to steradian.
    nbar        /=  (4. * np.pi / 41253.)
    nbar_wint   /=  (4. * np.pi / 41253.)

  elif band in ['u']: 
    ##  Malkan u-drops.
    mlim         =  24.5

    stats        =  usample_stats(mag=mlim)
    peakz        =  stats['Malkan']['z']

    alpha        =  stats['Malkan']['schechter']['alpha']
    Mstar        =  stats['Malkan']['schechter']['M_star']
    phi_star     =  stats['Malkan']['schechter']['phi_star']
    
    dzee         =  0.50             ##  0.61 / 2.
    peakz        =  stats['Malkan']['z']

    ##  Actual u-dropouts;  N per sq. deg. 
    nbar         =  projdensity(peakz - dzee / 2., peakz + dzee / 2., phi_star, Mstar, alpha, mlim=mlim, printit=True, completeness=None)
    nbar        /=  (4. * np.pi / 41253.)

    ##  With contamination. 
    nbar_wint    =  stats['Malkan']['nbar']    ##  per sq. deg.
    nbar_wint   /=  (4. * np.pi / 41253.)  ##  per steradian. 

  else:
    raise  ValueError('\n\nChosen sample is not available.\n\n')


  ##  Bias with z.
  drop_bz            =      get_dropoutbz(m=24.5) ## [24.5, 25.0, 25.5] 
  bz                 =  lambda z:  drop_bz(peakz)

  ##  LSST whitebook p(z).
  ##  pz             =  lambda z:  whitebook_pz(z, ilim = 25.30)
  
  ##  GOLDRUSH g-dropouts.
  ##  zee, pzee      =  get_gdropoutpz()
  
  ##  Get Rafelski based estimates on p(z).                                                                                                               
  midz, ps           =  Raf15_pz(droptype=band, field='UVUDF', depth='FULL', dz=0.1, no_lowz=True)
  pz                 =  interp1d(midz, ps, kind='linear', bounds_error=False, fill_value=0.0, assume_sorted=False)

  ##  Change in Ckg with p(z), b(z), nbar -> p'(z), b'(z) and nbar'.
  ##  We assume the z < 1 population of likely red galaxies has a bias of 2.04 at a mean z 0.87;
  bzz                =  lambda z:  bz(z)  if  z > 1.0  else  2.04 

  ##  pzz            =  lambda z:  whitebook_pz(z, ilim = 25.25)
  midzz, pss         =  Raf15_pz(droptype=band, field='UVUDF', depth='FULL', dz=0.1, no_lowz=False)
  pzz                =  interp1d(midzz, pss, kind='linear', bounds_error=False, fill_value=0.0, assume_sorted=False)

  ##  and the parameter Fisher matrix.                                                                                                                    
  fid_b1             =     bz(peakz)
  fid_sig8           =  sig8z(peakz)

  ##  Dry run.
  ##  cgg            =  Cgg(Pk_interps, Llls, zmin, zmax, pz, bz, zeff=True, bz2 = bz, survey_pz2 = pz)
  ##  ckg            =  Ckg(Pk_interps, Llls, zmin, zmax, pz, bz, zeff=True)
  ##  print(Llls, cgg, ckg)

  ##  {'gg': cgg + ngg, 'kg': ckg, 'gk': ckg, 'kk': ckk + nkk}
  Cls, Nls           =  get_allCls(Pk_interps, Llls, nbar, fsky, zmin, zmax, pz, bz, zeff=False, samplevar_lim=False)

  ##  Distorted Cls, i.e. with differing p(z), b(z) and nbar. 
  xCls, xNls         =  get_allCls(Pk_interps, Llls, nbar_wint, fsky, zmin, zmax, pzz, bzz, zeff=False, samplevar_lim=False)

  ##  Lll max cut;  Zel'dovich.                                                                                                                          
  LllMax             =  Lcutmax[round(peakz)][0]

  if plotit:
    pl.clf()

    latexify(columns=1, equal=True, fontsize=12)

    pl.axvline(LllMax, label='ZA limit.', c='k', alpha=0.5)

    pl.plot(Llls,       100. * np.abs(Cls['gg'] - xCls['gg']) / Cls['gg'],  label=r'$dC_{gg}/C_{gg} [\%]$')

    ##  kg is linear, gg is quad. in dNdz.
    pl.plot(Llls,       100. * np.abs(Cls['kg'] - xCls['kg']) / Cls['kg'],  label=r'$dC_{kg}/C_{kg} [\%]$')
    
    pl.xlim(50., 4.e3)
    pl.ylim(-1.,  30.)

    pl.xlabel(r'$L$')
    pl.legend(ncol=1, loc=1, frameon=False, handlelength=.6)

    plt.tight_layout()

    ##  pl.show()
    pl.savefig('plots/interloper_bias_cls_%s-drops.pdf' % band)  

  ##  The resulting parameter bias, (dsig8, db1).                                                                                                        
  nparam             =  2  
  dtheta             =  np.zeros(nparam)

  ##  dtheta stripped of iFab, eqn. (5.1) of https://arxiv.org/pdf/1706.03173.pdf 
  interim            =  np.zeros(nparam)

  ##  Parameter Fisher matrix.                                                                                                                        
  Fisher             =  np.zeros(nparam * nparam).reshape(nparam, nparam)
  iFisher            =  np.zeros_like(Fisher)

  ##  Per-L covariance;  e.g. 'ggkk' etc ... Sum of signal and noise.  eqn (14) of https://arxiv.org/pdf/1710.09465.pdf
  covs               =  get_ClsCov(Llls, Cls, Nls, fsky)

  ##  Covariance matrix of {kk, kg, gg} for given L = Llls[LL].
  nspectra           =  3
  Cov_Lll            =  np.zeros(nspectra * nspectra).reshape(nspectra, nspectra)

  ## and inverse.
  iCov_Lll           =  np.zeros_like(Cov_Lll) 
 
  ## Lll max cut;  Zel'dovich.                                                                                                                             
  Llls               =  Llls[Llls < LllMax]

  for LL, dummy in enumerate(Llls):
    dckg         =  Cls['kg'][LL] - xCls['kg'][LL]
    
    ##  Noise bias in Fisher expectation.
    dcgg         = (Cls['gg'][LL] + Nls['gg'][LL]) - (xCls['gg'][LL] + xNls['gg'][LL])

    ##  Construct D(LL) = [Ckk(LL), Ckg(LL), Cgg(LL)].  Noise bias?                                                                                  
    DL           =  np.array([Cls['kk'][LL], Cls['kg'][LL], Cls['gg'][LL]])

    for i, tracer in enumerate(['kk', 'kg', 'gg']):
      for j, ttracer in enumerate(['kk', 'kg', 'gg']):
        ##  Symmetric covariance between tracers at given LL, e.g. Cov(Ckk[LL], Ckg[LL]).
        Cov_Lll[i, j] = covs[tracer + ttracer][LL]

    ## .. and the inverse.  
    iCov_Lll = inv(Cov_Lll)

    for i, pp in enumerate(['s8', 'b1']):
      for j, ss in enumerate(['s8', 'b1']):
        ##  Per-L Fisher matrix.  Eqn. (17) of https://arxiv.org/pdf/1710.09465.pdf
        ##  Effectively working with samples in bins.  Fisher is a sum over modes.  Account for N>1 L modes per bin.
        Fisher[i, j] += nmodes[LL] * np.dot(dp_DL(pp, fid_sig8, fid_b1, DL), np.dot(iCov_Lll, dp_DL(ss, fid_sig8, fid_b1, DL)))

    ##  Eqn. (5.1) of https://arxiv.org/pdf/1706.03173.pdf;  No change to kk due to interloper population; Noise bias? 
    DDL     = np.array([0.0, dckg, dcgg])
    iCovDDL =   np.dot(iCov_Lll, DDL)
    
    for j, pparam in enumerate(['s8', 'b1']):
      ##  dmu_m / dp_beta * iCov(m,n) dDn.
      ##  Effectively working with samples in bins.  Fisher is a sum over modes.  Account for N>1 L modes per bin. 
      interim[j] += nmodes[LL] * np.dot(dp_DL(pparam, fid_sig8, fid_b1, DL), iCovDDL)
  
  ##  Invert Fisher. 
  iFisher = inv(Fisher)
  
  for i, param in enumerate(['s8', 'b1']):
    for j, pparam in enumerate(['s8', 'b1']):
      dtheta[i] += iFisher[i,j] * interim[j]
      
  biased_sig8 = fid_sig8 + dtheta[0]
  biased_b1   = fid_b1 + dtheta[1]

  print('\n\nFiducial sig8(z) and b1(z):  %.6lf +- %.6lf and %.6lf +- %.6lf.' % (fid_sig8, np.sqrt(iFisher[0,0]), fid_b1, np.sqrt(iFisher[1,1])))
  print('Biased   sig8(z) and b1(z):  %.6lf and %.6lf.' % (biased_sig8, biased_b1))
  print('Bias: %.6lf \t %.6lf [sigma]' % (np.abs(fid_sig8 - biased_sig8) / np.sqrt(iFisher[0,0]), np.abs(fid_b1 - biased_b1) / np.sqrt(iFisher[1,1])))

  ##  And plot contour ...                                                                                                                                 
  pl.clf()

  latexify(columns=1, equal=True, fontsize=12)

  fig    = plt.gcf()
  ax     = plt.gca()

  for mass_level, color, alpha in zip([0.99, 0.95, 0.68], ['b', 'b', 'b'], [0.2, 0.4, 0.6]):
    plot_ellipse(x_cent = fid_sig8, y_cent = fid_b1, ax = ax, cov = iFisher, mass_level = mass_level,\
                 fill=True, fill_kwargs={'alpha': alpha, 'c': color}, plot_kwargs={'c': color, 'alpha': 0.0})
  
  pl.plot(fid_sig8, fid_b1, 'w*', markersize=5, label=r'$(\sigma_8, b_1)$')
  pl.plot(biased_sig8, biased_b1, 'k*', markersize=5, alpha=0.4, label=r'$(\hat \sigma_8, \hat b_1)$')

  pl.legend()
  pl.xlabel(r'$\sigma_8(z=%.1lf)$' % peakz)
  pl.ylabel(r'$b_1(z=%.1lf)$'      % peakz)

  plt.tight_layout()

  ##  pl.show()
  pl.savefig('plots/interloper_bias_contours_%s-drops.pdf' % band)
  
  print("\n\nDone.\n\n")
