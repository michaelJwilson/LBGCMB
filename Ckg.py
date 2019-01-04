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


plt.style.use('ggplot')
mpl.rc('text', usetex = True)

params = get_params()

def Ckg(Pk_interps, Llls, zmin, zmax, survey_pz, bz, zeff = True):
  dz       = 0.01
  zs       = np.arange(zmin, zmax, dz)
  
  ##  Catch normalisation of p(z).  Added 03/01/19.                                                                                                                                                                                 
  ps       = survey_pz(zs)
  norm     = np.sum(ps) * dz
  ps      /= norm

  if zeff:
    ##  Calculate the mean redshift.
    zg      = np.sum(ps * zs) * dz / np.sum(dz * survey_pz(zs))

    chi_g   = comoving_distance(zg)

    k       = (Llls + 0.5) / chi_g

    ##  Limber and thin slice approximation.                                                                                                               
    result  = Pmm(Pk_interps, k, zg) * bz(zg) * lensing_kernel(zg) / chi_g**2.
    result /= params['h_100']**2.  ## Account for [h^-1 Mpc]^3 of Pmm and h^-1 Mpc of chi_g.                                                                

    return  result                 ## Dimensionless.                                                                                                        

  else:                            
     ## Assumes integral over dz slice.                                                                                       
     chis        = comoving_distance(zs)
     result      = np.zeros((len(zs), len(Llls)))

     for i, redshift in enumerate(zs):
       ks           = (Llls + 0.5) / chis[i]
       result[i,:]  = Pmm(Pk_interps, ks, redshift) * bz(redshift)

     prefactor   =  lensing_kernel(zs) * sliced_pz(zs, zmin, zmax, survey_pz)/chis**2.
     integrand   =  prefactor[:, None] * result

     integrand  /=  params['h_100'] **2.  ## Account for [h^-1 Mpc]^3 of Pmm and [h^-1 Mpc] of chi.                                                         

     ## Integrated over z, Ckg(L).                                                                                                                        
     result      = simps(integrand, dx = zs[1] - zs[0], axis=0)

     return  result

def var_Ckg(Pk_interps, lensCl_interps, nolensCl_interps, Llls, zmin, zmax, survey_pz, bz, nbar, fsky,\
            ngg=None, nkk=None, thetab=1., DeltaT=1., zeff = True, modesample_lim=False, iterative=False, pickleit=True):

  ckk = Ckk(Pk_interps, Llls)
  cgg = Cgg(Pk_interps, Llls, zmin, zmax, survey_pz, bz, zeff = zeff)
  ckg = Ckg(Pk_interps, Llls, zmin, zmax, survey_pz, bz, zeff = zeff)

  if  modesample_lim:
    ngg   = np.zeros_like(Llls)

  else:
    if ngg is None:
      ngg = Ngg(Llls, zmin, zmax, survey_pz, nbar)

  if nkk is None:
    ##  Ordering within terms is important. 
    nkk   = Nkk(lensCl_interps, nolensCl_interps, Llls, ['TT', 'TE', 'EE', 'EB'], thetab, DeltaT, iterative=iterative, pickleit=pickleit)

  num     = (ckk + nkk) * (cgg + ngg) + ckg**2.   ## Gaussian approximation.                                                                               
  den     = (2. * Llls + 1.) * fsky

  return  num / den

  
if __name__ == "__main__":
  import  pylab              as      pl

  from    prep_Llls          import  prep_Llls
  from    Gaussian_pz        import  Gaussian_pz
  from    snr                import  snr
  from    pmh                import  get_PkInterps, linz_bz
  from    prep_camb          import  CAMB, Clxy
  from    bolometers         import  bolometers
  from    completeness       import  get_dropoutpz
  # from  dropouts.reddy     import  samplestats
  from    specs              import  samplestats
  from    ilim               import  get_nbar_nocontam
  from    schmittfull        import  ss_pz
  from    scipy.interpolate  import  interp1d  
  from    bz                 import  get_dropoutbz


  print("\n\nWelcome to Ckg.\n\n")

  cmbexp                          =  'Planck'
  fsky, thetab, DeltaT, iterative =  bolometers[cmbexp]['fsky'],   bolometers[cmbexp]['thetab'],\
                                     bolometers[cmbexp]['DeltaT'], bolometers[cmbexp]['iterative']

  ## Prepare pycamb module; linear, non-linear matter P(k) and Cls.                                                                                        
  cambx                =  CAMB()

  NLlls, Llls, nmodes  =  prep_Llls(NLlls = 60, Lmin = 50., Lmax = 5000., log10=True)
  Pk_interps           =  get_PkInterps(cambx)

  ## No Detector noise -- this is handled by Clxy of prep_camb.                                                                                            
  (lensCl_interps, nolensCl_interps) = cambx.get_Cls()

  ## Galaxies per sq. degree.                                                                                                                              
  band                 =  'g'
  stats                =  samplestats()
  stats                =  get_nbar_nocontam(band, depth='W', printit=False)

  nbar                 =  stats[band]['nbar_nointerlopers']
  peakz                =  stats[band]['z']
  bz                   =  linz_bz ## get_dropoutbz()
  
  zmin                 =  peakz - 2.00
  zmax                 =  peakz + 2.00

  ##  Effectively overwrites hard z limits above.
  zee, pzee            =  get_dropoutpz()
  pz                   =  interp1d(zee, pzee, kind='linear', bounds_error=False, fill_value=0.0, assume_sorted=False)

  ## Schmittfull and Seljak (2017).
  ## pz, nbar          =  ss_pz()

  ckg                 =  Ckg(Pk_interps, Llls, zmin, zmax, pz, bz, zeff=True)
  nkk                 =  Nkk(lensCl_interps, nolensCl_interps, Llls, terms=['TT', 'TE', 'EE', 'EB'],\
                             thetab=thetab, DeltaT=DeltaT, iterative=iterative, pickleit=False)

  vkg                 =  var_Ckg(Pk_interps, lensCl_interps, nolensCl_interps, Llls, zmin, zmax, pz, bz, nbar, fsky,\
                                 nkk=nkk, iterative=iterative, pickleit=True)
  
  print('\n\nTotal S/N:  %.1lf' % snr(Llls, ckg, vkg, nmodes, lmax=1.e5))

  ax                  = pl.gca()
  ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))

  pl.loglog(Llls, ckg, c='g')
  pl.errorbar(Llls, ckg, np.sqrt(vkg), c='g')

  pl.xlim(4.e1,   5.e3)
  pl.ylim(1.e-9, 1.e-4)

  pl.xlabel(r'$L$')
  pl.ylabel(r'$C_{\kappa g}(L)$')

  pl.savefig('plots/Ckg.pdf', bbox_inches='tight')
  
  print("\n\nDone.\n\n")
