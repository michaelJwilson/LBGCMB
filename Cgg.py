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
from    pmh                import  Pmm, get_PkInterps, linz_bz
from    reddy              import  samplestats as reddy_samplestats
from    sliced_pz          import  sliced_pz
from    bz                 import  get_dropoutbz
from    utils              import  latexify
from    nbar               import  projdensity
from    zeldovich_Lmax     import  Lcutmax
from    matplotlib.patches import  Rectangle


latexify(fig_width=None, fig_height=None, columns=1, equal=True, fontsize=10)

params = get_params()
         
def Cgg(Pk_interps, Llls, zmin, zmax, survey_pz, bz, survey_pz2=None, bz2=None, zeff=True):
  dz          = 0.001
  zs          = np.arange(zmin, zmax + dz, dz)

  ##  Catch normalisation of p(z).  Added 03/01/19.
  ps          = survey_pz(zs)
  norm        = np.sum(ps) * dz
  ps         /= norm

  chis        = comoving_distance(zs)

  if bz2 is None:
     bz2      = bz

  if survey_pz2 is None:
     survey_pz2 = survey_pz
     ps2        = np.copy(ps)

  else:
    ps2         = survey_pz2(zs)
    norm        = np.sum(ps2) * dz
    ps2        /= norm

  if zeff: 
    ''' 
    Compute Cgg in the z_eff and first-order Limber approximations 
    for a slice of galaxies at zz of width dz;  eqn. (5) of 
    https://arxiv.org/pdf/1511.04457.pdf
    '''
    zz        =  np.sum(dz * ps * zs) / (np.sum(ps) * dz)                      ##  Normalisation should already be unity.  
    ks        = (Llls + 0.5) / comoving_distance(zz)                           ##  For the Phh evaluation in the integral, we take a zeff approx. 
                                                                               ##  i.e. \int dz .... Phh(zeff).

    result    =  Pmm(Pk_interps, ks, zz) * bz(zz) * bz2(zz)
    result    =  np.broadcast_to(result, (len(zs), len(Llls)))                 ##  Broadcast to each redshift.  

  else: 
    result    = np.zeros((len(zs), len(Llls)))

    for i, z in enumerate(zs):
      ks          = (Llls + 0.5) / chis[i]                                     ##  Evaluate Phh for each z and k(ell, z).
      result[i,:] =  Pmm(Pk_interps, ks, z) * bz(z) * bz2(z)

  ##  Accounts for both spatial to angular mapping as function of z;  i.e. k = (l. + 0.5) / chi(z) and redshift evolution of P_mm(k).
  prefactor   = (cosmo.H(zs).value / const.c.to('km/s').value) * sliced_pz(zs, zmin, zmax, survey_pz) * sliced_pz(zs, zmin, zmax, survey_pz2) / chis**2.
  
  integrand   =  prefactor[:, None] * result
  integrand  /=  params['h_100']

  return  simps(integrand, dx = dz, axis=0)                                    ## integral over z. 

def var_Cgg(Llls, zmin, zmax, survey_pz, bz, nbar, fsky, cgg = None, Pk_interps=None, samplevar_lim=False):
  if cgg is None:    
    cgg  = Cgg(Pk_interps, Llls, zmin, zmax, survey_pz, bz, zeff=True)

  if samplevar_lim:
    num  = cgg

  else:
    num  = cgg + Ngg(Llls, zmin, zmax, survey_pz, nbar)

  num    =  2. * num ** 2.  
  den    = (2. * Llls + 1.) * fsky

  return  num / den

def Ngg(Llls, zmin, zmax, survey_pz, nbar):  
  '''                                                                                                                                            
  Poisson noise on Cgg equal to (1. / angular number density).                                                                   
  '''   
  from  pz2nbar  import  pz_slice, nbar_convert


  spz   = pz_slice(zmin, zmax, survey_pz)

  ##  Assumes input nbar is 1 / deg2 (the default unit) Converts to  per steradian. 
  nbar  = nbar_convert(nbar, unit='str')

  return  np.ones_like(Llls) / (nbar * spz)

def maglim_ax(Llls, cgg, ax, band = 'g', decband='i'):
  from  pz2nbar  import  nbar_convert


  root       = os.environ['LBGCMB']

  data       = np.loadtxt(root + "/dropouts/schechter/dat/schechter_estimate_%s_dropouts.txt" % band)
  ms         = data[:,0]
  Ns         = data[:,1]

  Ns         = nbar_convert(Ns, unit='str')

  ms_interp  = interp1d(Ns, ms,      kind='linear', bounds_error=True, fill_value=0.0, assume_sorted=False) 
  iNs_interp = interp1d(ms, 1. / Ns, kind='linear', bounds_error=True, fill_value=0.0, assume_sorted=False)

  ##  Add second mag_axis.                                                                                                                               
  ax2        = ax.twinx()
  ymin, ymax = ax.get_ylim()

  ax2.loglog(Llls, cgg, alpha=0.0)

  ax2.set_ylim([ymin, ymax])

  ## ticks   = [5.e-8, 1.e-7, 5.e-7, 10.e-7, 50.e-7, 100.e-7, 500.e-7, 1000.e-7]
  ## strings = ['%.1lf' % ms_interp(1. / x) for x in ticks]

  tms        =  np.arange(23.5, 26.5, 0.5)

  strings    = ['%.1lf' % x   for x in tms]
  tms        = [iNs_interp(x) for x in tms]

  print(strings)
  print(tms)

  ax2.set_yticks(tms)
  ax2.set_yticklabels(strings)

  ax2.set_ylabel(r'$%s_{{\rm{AB}}}$' % decband)

  ax2.tick_params(axis='y', which='minor', color='w')

  return  ax2


if __name__ == "__main__":
  import  pylab              as      pl

  from    prep_Llls          import  prep_Llls
  from    Gaussian_pz        import  Gaussian_pz
  from    prep_camb          import  CAMB
  ##  from    schmittfull        import  ss_pz
  from    snr                import  snr
  from    completeness       import  get_dropoutpz       as  get_gdropoutpz
  from    Malkan.specs       import  samplestats         as  usample_stats                                                                                 
  from    specs              import  samplestats         as  gsample_stats
  from    ilim               import  get_nbar_nocontam
  from    scipy.interpolate  import  interp1d
  from    bolometers         import  bolometers
  from    lensing            import  Ckk, var_Ckk
  from    Nkk                import  Nkk
  from    Ckg                import  Ckg, var_Ckg 


  print("\n\nWelcome to Cgg.\n\n")

  ##  Prepare pycamb module; linear, non-linear matter P(k) and Cls.                                                                                     
  cambx                              =  CAMB()
  Pk_interps                         =  get_PkInterps(cambx)

  ##  No Detector noise -- this should be handled by Clxy of prep_camb.                                                                                  
  (lensCl_interps, nolensCl_interps) =  cambx.get_Cls()

  NLlls, Llls, nmodes                =  prep_Llls(NLlls = 60, Lmin = 50., Lmax = 5000., log10=True)

  cmbexp                             =  'CMBS4'

  fsky, thetab, DeltaT, iterative    =  bolometers[cmbexp]['fsky'],   bolometers[cmbexp]['thetab'],\
                                        bolometers[cmbexp]['DeltaT'], bolometers[cmbexp]['iterative']
  
  ##  band      =  'Malkan'              ## Reddy u-drops.
  band  =  'g'
  
  
  ##  Set (no interloper) nbar, b(z) and p(z).                                                                                                         
  
  if  band == 'g':
    stats        =  gsample_stats()
    stats        =  get_nbar_nocontam(band, depth='W', printit=False)
  
    ##  Effectively overwrites hard z limits above.                                                                                                 
    zee, pzee    =  get_gdropoutpz()
    pz           =  interp1d(zee, pzee, kind='linear', bounds_error=False, fill_value=0.0, assume_sorted=False)

    nbar         =  stats[band]['nbar_noint']
    peakz        =  stats[band]['z']

    decband      =  'i'
    colors       =  ['darkgreen', 'limegreen', 'g']
  
  elif band ==  'Malkan': 
    ##  Malkan u-drops.
    stats        =  usample_stats()

    peakz        =  stats[band]['z']

    alpha        =  stats[band]['schechter']['alpha']
    Mstar        =  stats[band]['schechter']['M_star']
    phi_star     =  stats[band]['schechter']['phi_star']
    
    dzee         =  0.61 / 2.
    peakz        =  stats[band]['z']
    nbar         =  projdensity(peakz - dzee / 2., peakz + dzee / 2., phi_star, Mstar, alpha, mlim=24.5, printit = True, completeness=None)
  
    pz           =  Gaussian_pz

    decband      =  'R'
    colors       =  ['darkblue', 'deepskyblue', 'b']

  else:
    raise ValueError('\n\nChosen band is not available.\n\n')

  ## Galaxies per sq. degree.                                                                                                                           
  ## pz, nbar  =  ss_pz()  

  ##  Bias with z.
  bz           =  linz_bz
  ##  bz       =  get_dropoutbz()

  ##  hard p(z) limits.
  zmin         =  peakz - 2.00
  zmax         =  peakz + 2.00

  ##
  pl.clf() 

  ##  Cgg.
  cgg          =      Cgg(Pk_interps, Llls, zmin, zmax, pz, bz, zeff=True, bz2 = bz, survey_pz2 = pz)
  ngg          =      Ngg(Llls, zmin, zmax, pz, nbar)
  vgg          =  var_Cgg(Llls, zmin, zmax, pz, bz, nbar, fsky, samplevar_lim=False, cgg = cgg)

  print('\n\nTotal Cgg S/N:  %.1lf' % snr(Llls, cgg, vgg, nmodes, lmax=1.e5))

  pl.loglog(Llls, cgg, colors[0], label=r'$C_{gg}$')

  pl.axhline(y = ngg[0], xmin = 0., xmax = 1.e4, c=colors[0], alpha=0.4, label=r'$N_{gg}$') 

  ax           = pl.gca()                                                                                                                            
  #ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))

  pl.errorbar(Llls, cgg, np.sqrt(vgg), c=colors[0])
  
  ## Ckg.
  ckg                 =  Ckg(Pk_interps, Llls, zmin, zmax, pz, bz, zeff=True)
  nkk                 =  Nkk(lensCl_interps, nolensCl_interps, Llls, terms=['TT', 'TE', 'EE', 'EB'],\
                             thetab=thetab, DeltaT=DeltaT, iterative=iterative, pickleit=False)

  vkg                 =  var_Ckg(Pk_interps, lensCl_interps, nolensCl_interps, Llls, zmin, zmax, pz, bz, nbar, fsky,\
                                 nkk=nkk, iterative=iterative, pickleit=True)

  pl.loglog(Llls,   ckg, c=colors[1], label=r'$C_{\kappa g}$')
  pl.errorbar(Llls, ckg, np.sqrt(vkg), c=colors[1])
  
  ## Ckk
  ckk                  = Ckk(Pk_interps, Llls)
  nkk                  = Nkk(lensCl_interps, nolensCl_interps, Llls, terms=['TT', 'TE', 'EE', 'EB'], thetab=thetab, DeltaT=DeltaT, iterative=iterative,\
                             pickleit=True)

  vkk                  = var_Ckk(Llls, fsky, nkk, Pk_interps, samplevar_lim=False)

  pl.loglog(Llls, ckk,   'k', label=r'$C_{\kappa \kappa}$')
  pl.loglog(Llls, nkk,   'k', alpha=0.4, label=r'$N_{\kappa \kappa}$')

  pl.errorbar(Llls, ckk, np.sqrt(vkk), c='k')

  ##  Plot Lmax modelling limit for this redshift. 
  ax                   =  pl.gca()
  Lmax                 =  Lcutmax[np.round(peakz)][0]
  

  ax.add_patch(Rectangle((Lmax, 0.), 5.e3,    1.e6, color=colors[2], alpha=0.3))

  ##  Plot SO noise curve. 
  for cmbexp in ['SO', 'Planck']:
    fsky, thetab, DeltaT, iterative    =  bolometers[cmbexp]['fsky'],   bolometers[cmbexp]['thetab'],\
                                          bolometers[cmbexp]['DeltaT'], bolometers[cmbexp]['iterative']

    nkk                                =  Nkk(lensCl_interps, nolensCl_interps, Llls, terms=['TT', 'TE', 'EE', 'EB'],\
                                          thetab=thetab, DeltaT=DeltaT, iterative=iterative, pickleit=False)

    pl.loglog(Llls, nkk, 'k', alpha=0.4, label=r'')
  
  pl.xlim(5.e1,   5.e3)
  pl.ylim(1.e-9, 3.e-5)

  pl.xlabel(r'$L$')
  
  pl.legend(loc = 2, ncol=3, frameon=False, facecolor='w')

  ## Mag axis.                                                                                                            
  ax2 = maglim_ax(Llls, cgg, ax, band = band, decband=decband)
  ## ax2.set_axis_on()

  ax2.spines['bottom'].set_color('black')
  ax2.spines['top'].set_color('black')
  ax2.spines['left'].set_color('black')
  ax2.spines['right'].set_color('black')
  ax2.grid(False)

  pl.savefig('plots/%sCgg.pdf' % band, bbox_inches='tight')
  
  print("\n\nDone.\n\n")
