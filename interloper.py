import  os
import  json
import  pickle
import  matplotlib          as      mpl
import  matplotlib.pyplot   as      plt
import  numpy               as      np
import  astropy.units       as      u
  
from    params              import  get_params
from    pmh                 import  get_PkInterps
from    utils               import  latexify
from    nbar                import  projdensity
from    zeldovich_Lmax      import  Lcutmax
from    Cgg                 import  Cgg, var_Cgg, Ngg
from    numpy.linalg        import  inv
from    fisher_contour      import  plot_ellipse
from    planck18_bao        import  get_sig8z
from    interlopers.pz      import  get_pz             as gp_pz


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
  from    prep_camb          import  CAMB
  from    scipy.interpolate  import  interp1d
  from    bolometers         import  bolometers
  from    lensing            import  Ckk, var_Ckk
  from    Nkk                import  Nkk
  from    Ckg                import  Ckg, var_Ckg 
  from    schechter.gen_pz   import  peakz               as  _peakz
  from    schechter.get_shot import  get_shot
  from    schechter.get_pz   import  get_pz
  from    get_bz             import  bz_callmodel


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
  

  ##  
  band       =   'u'

  setup      = {'BX': {'colors': ['goldenrod', 'tan',         'y'], 'maglim': 25.5, 'decband': 'R'},\
                 'u': {'colors': ['darkblue',  'deepskyblue', 'b'], 'maglim': 24.6, 'decband': 'R'},\
                 'g': {'colors': ['darkgreen', 'limegreen',   'g'], 'maglim': 25.8, 'decband': 'i'},\
                 'r': {'colors': ['darkred',   'indianred',   'r'], 'maglim': 25.8, 'decband': 'z'}}

  mlim       =  setup[band]['maglim']

  ##  pz     =  get_pz(band)

  ##  Interloper GP rewrite. 
  _, _, _, _, ngp = gp_pz(band, 'Full', 5.5)
  pz              = np.vectorize(ngp)
  
  ## 
  bz         =  lambda z:  bz_callmodel(z, mlim)

  nbar       =  get_shot(band, mlim)
  nbar      /= (4. * np.pi / 41253.)  ##  Sq. deg. to steradian. 

  ## 
  int_frac   =  0.06
  nbar_wint  =  (1.0 + int_frac) * nbar 

  decband    =  setup[band]['decband']
  colors     =  setup[band]['colors']

  peakz      =  _peakz(pz)

  ##  Hard p(z) limits. 
  zmin       =   0.01
  zmax       =  10.00
  
  ##  LSST whitebook p(z).
  ##  pz             =  lambda z:  whitebook_pz(z, ilim = 25.30)
  
  ##  Change in Ckg with p(z), b(z), nbar -> p'(z), b'(z) and nbar'.
  ##  We assume the z < 1 population of likely red galaxies has a bias of 2.04 at a mean z 0.87;
  ##  bzz            =  lambda z:  bz(z)  if  z > 1.0  else  2.04 
  bzz                =  lambda z:  bz(z)  if  z > 3.0  else (1. + z)

  ##  Interloper GP rewrite.                                                                                                                                                    
  _, _, _, _, _ngp   = gp_pz(band, 'Degraded', 5.5)
  pzz                = np.vectorize(_ngp)

  ##  and the parameter Fisher matrix.                                                                                                                    
  fid_b1             =     bz(peakz)
  fid_sig8           =  sig8z(peakz)

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

    pl.plot(Llls,  100. * np.abs(Cls['gg'] - xCls['gg']) / Cls['gg'],  label=r'$|dC_{gg}/C_{gg}| [\%]$')

    ##  kg is linear, gg is quad. in dNdz.
    pl.plot(Llls,  100. * np.abs(Cls['kg'] - xCls['kg']) / Cls['kg'],  label=r'$|dC_{kg}/C_{kg}| [\%]$')
    
    pl.xlim(50., 2.5e3)
    ##  pl.ylim(100., 450.)

    pl.xlabel(r'$L$')
    pl.legend(ncol=1, loc=2, frameon=False, handlelength=.6)

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

  ##  and inverse.
  iCov_Lll           =  np.zeros_like(Cov_Lll) 
 
  ##  Lll max cut;  Zel'dovich.                                                                                                                             
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
        ##  Effectively working with samples in bins.  Fisher is a sum over modes.  Account for N > 1 L modes per bin.
        Fisher[i, j] += nmodes[LL] * np.dot(dp_DL(pp, fid_sig8, fid_b1, DL), np.dot(iCov_Lll, dp_DL(ss, fid_sig8, fid_b1, DL)))

    ##  Eqn. (5.1) of https://arxiv.org/pdf/1706.03173.pdf;  No change to kk due to interloper population. 
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

  ##  
  biased_sig8 = fid_sig8 + dtheta[0]
  biased_b1   = fid_b1   + dtheta[1]

  print('\n\nFiducial sig8(z) and b1(z):  %.6lf +- %.6lf and %.6lf +- %.6lf.' % (fid_sig8, np.sqrt(iFisher[0,0]), fid_b1, np.sqrt(iFisher[1,1])))
  print('Biased   sig8(z) and b1(z):  %.6lf and %.6lf.' % (biased_sig8, biased_b1))
  print('Bias: %.6lf \t %.6lf [sigma]' % (np.abs(fid_sig8 - biased_sig8) / np.sqrt(iFisher[0,0]), np.abs(fid_b1 - biased_b1) / np.sqrt(iFisher[1,1])))

  ##  Save results as json. 
  output = {'peakz': peakz, 'fid_sig8': np.float(fid_sig8), 'fid_b1': np.float(fid_b1), 'biased_sig8': biased_sig8, 'biased_b1': biased_b1, 'iFisher': iFisher.tolist(),\
            'intfrac': int_frac}

  with open('dat/result4interloper_%s_intfrac_%.2lf.json' % (band, int_frac), 'w') as fp:
    json.dump(output, fp, sort_keys=True, indent=4)

  print("\n\nDone.\n\n")
