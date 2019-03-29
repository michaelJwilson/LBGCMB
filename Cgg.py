import  os
import  numpy              as      np
import  astropy.units      as      u

from    cosmo              import  cosmo
from    params             import  get_params
from    utils              import  comoving_distance
from    astropy            import  constants as const
from    scipy.integrate    import  simps
from    pmh                import  Pmm, get_PkInterps, linz_bz
from    sliced_pz          import  sliced_pz
from    pz2nbar            import  pz_slice, nbar_convert


params = get_params()
         
def Cgg(Pk_interps, Llls, zmin, zmax, survey_pz, bz, survey_pz2=None, bz2=None, zeff=True):
  dz            = 0.01
  zs            = np.arange(zmin, zmax + dz, dz)

  ##  Catch normalisation of p(z).  Added 03/01/19.
  ps            = survey_pz(zs)
  norm          = np.sum(ps) * dz
  ps           /= norm

  chis          = comoving_distance(zs)

  if bz2 is None:
     bz2        = bz

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
    zz        =  dz * np.sum(ps * zs) / (np.sum(ps) * dz)                      ##  Normalisation should already be unity.  
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
  prefactor   = (cosmo.H(zs).value / const.c.to('km/s').value) * sliced_pz(zs, zmin, zmax, survey_pz)\
                                                               * sliced_pz(zs, zmin, zmax, survey_pz2)\
                                                               / chis**2.
  
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

  spz   = pz_slice(zmin, zmax, survey_pz)

  ##  Assumes input nbar is per sq. deg. (the default unit) to convert to per steradian. 
  nbar  = nbar_convert(nbar, unit='str')

  return  np.ones_like(Llls) / (nbar * spz)

def maglim_ax(Llls, cgg, ax, band = 'g', decband='i'):
  from  pz2nbar  import  nbar_convert

  ms, Ns      = get_shot(type=band, mlim=None)
  Ns          = nbar_convert(Ns, unit='str')

  iNs_interp  = interp1d(ms, 1. / Ns, kind='linear', bounds_error=True, assume_sorted=False)

  ##  Add second mag_axis.                                                                                                                        
  ax2         = ax.twinx()
  ymin, ymax  = ax.get_ylim()

  ax2.loglog(Llls, cgg, alpha=0.0)
    
  tms         =  np.arange(23.0, 26.0, 0.25)

  strings     = ['%.1lf' % x   for x in tms]
  tms         = [iNs_interp(x) for x in tms]

  ax2.set_yticks(tms)
  ax2.set_yticklabels(strings)
  
  ax2.set_ylabel(r'$%s_{{\rm{AB}}}$' % decband)
  ax2.tick_params(axis='y', which='minor', color='w')
  
  ax2.set_ylim([ymin, ymax])

  return  ax2


if __name__ == "__main__":
  import  pylab                    as      pl
  import  matplotlib               as      mpl
  import  matplotlib.pyplot        as      plt
  
  from    prep_Llls                import  prep_Llls
  from    prep_camb                import  CAMB
  from    scipy.interpolate        import  interp1d
  from    bolometers               import  bolometers
  from    lensing                  import  Ckk, var_Ckk
  from    Nkk                      import  Nkk
  from    Ckg                      import  Ckg, var_Ckg 
  from    utils                    import  latexify
  from    zeldovich_Lmax           import  Lcutmax
  from    matplotlib.patches       import  Rectangle
  from    schechter.gen_pz         import  peakz            as  _peakz
  from    schechter.get_shot       import  get_shot
  from    schechter.get_pz         import  get_pz
  from    get_bz                   import  bz_callmodel


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
  
  band       =  'BX'                      

  setup      = {'BX': {'colors': ['goldenrod', 'tan',         'y'], 'maglim': 25.5, 'decband': 'R'},\
                 'u': {'colors': ['darkblue',  'deepskyblue', 'b'], 'maglim': 24.6, 'decband': 'i'},\
                 'g': {'colors': ['darkgreen', 'limegreen',   'g'], 'maglim': 25.8, 'decband': 'i'},\
                 'r': {'colors': ['darkred',   'indianred',   'r'], 'maglim': 25.8, 'decband': 'z'}}
  
  mlim       =  setup[band]['maglim']

  pz         =  get_pz(band)
  bz         =  lambda z:  bz_callmodel(z, mlim)

  nbar       =  get_shot(band, mlim)              

  decband    =  setup[band]['decband']
  colors     =  setup[band]['colors']

  peakz      =  _peakz(pz)

  ##  Hard p(z) limits.  
  zmin       =  peakz - 2.00
  zmax       =  peakz + 2.00

  ##
  pl.clf() 

  latexify(fig_width=None, fig_height=None, columns=1, equal=True, fontsize=10)
  
  ##  Cgg.
  cgg          =      Cgg(Pk_interps, Llls, zmin, zmax, pz, bz, zeff=False, bz2 = bz, survey_pz2 = pz)
  ngg          =      Ngg(Llls, zmin, zmax, pz, nbar)
  vgg          =  var_Cgg(Llls, zmin, zmax, pz, bz, nbar, fsky, samplevar_lim=False, cgg = cgg)
 
  print(nbar, ngg[0], zmin, zmax, pz_slice(zmin, zmax, pz))
 
  pl.loglog(Llls, cgg, colors[0], label=r'$C_{gg}$')
  pl.axhline(y = ngg[0], xmin = 0., xmax = 1.e4, c=colors[0], alpha=0.4, label=r'$N_{gg}$') 

  ax           = pl.gca()                                                                                                                            
  #ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))

  pl.errorbar(Llls, cgg, np.sqrt(vgg), c=colors[0])
  
  ##  Ckg.
  ckg                 =  Ckg(Pk_interps, Llls, zmin, zmax, pz, bz, zeff=True)
  nkk                 =  Nkk(lensCl_interps, nolensCl_interps, Llls, terms=['TT', 'TE', 'EE', 'EB'],\
                             thetab=thetab, DeltaT=DeltaT, iterative=iterative, pickleit=False)

  vkg                 =  var_Ckg(Pk_interps, lensCl_interps, nolensCl_interps, Llls, zmin, zmax, pz, bz, nbar, fsky,\
                                 nkk=nkk, iterative=iterative, pickleit=True)

  pl.loglog(Llls,   ckg, c=colors[1], label=r'$C_{\kappa g}$')
  pl.errorbar(Llls, ckg, np.sqrt(vkg), c=colors[1])

  ##  Plot Lmax modelling limit for this redshift.                                                                                      
  ax                   =  pl.gca()
  Lmax                 =  Lcutmax[np.round(peakz)][0]

  ax.add_patch(Rectangle((Lmax, 0.), 5.e3,    1.e6, color=colors[2], alpha=0.3))
  
  ##  Ckk
  ckk                  = Ckk(Pk_interps, Llls)
  nkk                  = Nkk(lensCl_interps, nolensCl_interps, Llls, terms=['TT', 'TE', 'EE', 'EB'], thetab=thetab, DeltaT=DeltaT,\
                             iterative=iterative, pickleit=True)

  vkk                  = var_Ckk(Llls, fsky, nkk, Pk_interps, samplevar_lim=False)

  pl.loglog(Llls, ckk,   'k', label=r'$C_{\kappa \kappa}$')
  pl.errorbar(Llls, ckk, np.sqrt(vkk), c='k')

  pl.loglog(Llls[Llls < Lmax], nkk[Llls < Lmax],   'k', alpha=0.4, label=r'$N_{\kappa \kappa}$')

  ##  Plot SO and Planck noise curves. 
  for cmbexp in ['SO', 'Planck']:
    fsky, thetab, DeltaT, iterative    =  bolometers[cmbexp]['fsky'],   bolometers[cmbexp]['thetab'],\
                                          bolometers[cmbexp]['DeltaT'], bolometers[cmbexp]['iterative']

    nkk                                =  Nkk(lensCl_interps, nolensCl_interps, Llls, terms=['TT', 'TE', 'EE', 'EB'],\
                                          thetab=thetab, DeltaT=DeltaT, iterative=iterative, pickleit=False)

    pl.loglog(Llls[Llls < Lmax], nkk[Llls < Lmax], 'k', alpha=0.4, label=r'')
  
  pl.xlim(5.e1,   5.e3)
  pl.ylim(1.e-9, 3.e-5)

  pl.xlabel(r'$L$')
  
  pl.legend(loc = 2, ncol=3, frameon=False, handlelength=1.5, columnspacing=1.5)

  ##  Mag axis.                                                                                                            
  ax2 = maglim_ax(Llls, cgg, ax, band = band, decband=decband)

  ax2.spines['bottom'].set_color('black')
  ax2.spines['top'].set_color('black')
  ax2.spines['left'].set_color('black')
  ax2.spines['right'].set_color('black')
  ax2.grid(False)

  plt.tight_layout()

  ##  pl.show()
  pl.savefig('plots/%sCgg.pdf' % band, bbox_inches='tight')
  
  print("\n\nDone.\n\n")
