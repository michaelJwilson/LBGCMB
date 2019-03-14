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
from    sliced_pz          import  sliced_pz
from    utils              import  latexify
from    scipy.interpolate  import  interp1d
from    schmittfull_nz     import  get_ss17_samples


params   = get_params()

def snr(Llls, Cls, Vars, nmodes, lmax=1.e5):
  result  = (Cls[Llls <= lmax] / np.sqrt(Vars[Llls <= lmax])) ** 2.
  result  = np.sum(nmodes[Llls <= lmax] * result)

  return  np.sqrt(result)

def zeros(arg):
  return  np.zeros_like(arg)


if __name__ == '__main__':
  import  pylab               as      pl

  from    prep_Llls           import  prep_Llls
  from    pmh                 import  Pmm, get_PkInterps, linz_bz
  from    prep_camb           import  CAMB
  from    scipy.interpolate   import  interp1d
  from    prep_camb           import  CAMB
  from    bolometers          import  bolometers
  from    lensing             import  Ckk, var_Ckk
  from    Nkk                 import  Nkk
  from    Ckg                 import  Ckg, var_Ckg
  from    bolometers          import  bolometers
  from    zeldovich_Lmax      import  Lcutmax
  from    pylab               import  rcParams
  from    schechter.gen_pz    import  peakz            as  _peakz
  from    schechter.get_shot  import  get_shot
  from    schechter.get_pz    import  get_pz
  from    get_bz              import  bz_callmodel


  print('\n\nWelcome to snr.\n\n')
  
  band       =   'g'

  setup      = {'BX': {'colors': ['goldenrod', 'tan',         'y'], 'maglim': 25.5, 'decband': 'R'},\
                 'u': {'colors': ['darkblue',  'deepskyblue', 'b'], 'maglim': 25.5, 'decband': 'R'},\
                 'g': {'colors': ['darkgreen', 'limegreen',   'g'], 'maglim': 25.5, 'decband': 'i'},\
                 'r': {'colors': ['darkred',   'indianred',   'r'], 'maglim': 25.5, 'decband': 'z'}}

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


  Lmax                 =  Lcutmax[np.round(peakz)][0] 

  ##  Prepare pycamb module; linear, non-linear matter P(k) and Cls.                                                                        
  cambx                =  CAMB()
  NLlls, Llls, nmodes  =  prep_Llls(NLlls = 60, Lmin = 50., Lmax = 5000., log10=True)

  Pk_interps           =  get_PkInterps(cambx)

  ##  No Detector noise -- this is handled by Clxy of prep_camb.                                                                                   
  (lensCl_interps, nolensCl_interps) = cambx.get_Cls()

  ##  Get low-z delensing efficiency. 
  ns, ps, bs, ss       =  get_ss17_samples(nolsst=True)  ##  <\bar n>, p(z), b(z).

  rho                  =  np.loadtxt('rho/' + '_'.join(s for s in ss) + '.txt')
  rho                  =  interp1d(rho[:,0], rho[:,1], kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)

  ##  Dilution factors.
  dfactors             =  np.array([1.e-5, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 1.e1, 5.e1, 3.e3])
  colors               =  plt.rcParams['axes.prop_cycle'].by_key()['color']

  ##  
  latexify(fig_width=None, fig_height=None, columns=1, equal=True, fontsize=10)
  
  for ii, cmbexp in enumerate(bolometers): 
    ckg                   =  Ckg(Pk_interps, Llls, zmin, zmax, pz, bz, zeff=True)

    fsky, thetab, DeltaT, iterative = bolometers[cmbexp]['fsky'],   bolometers[cmbexp]['thetab'],\
                                      bolometers[cmbexp]['DeltaT'], bolometers[cmbexp]['iterative']

    nkk                   =  Nkk(lensCl_interps, nolensCl_interps, Llls, terms=['TT', 'TE', 'EE', 'EB'],\
                                 thetab=thetab,  DeltaT=DeltaT,    iterative=iterative, pickleit=False)

    for rho, ls, label in zip([zeros, rho], ['-', '--'], [cmbexp, '']):
      results             =  []

      for dfactor in dfactors:    
        dnbar             = dfactor * nbar
 
        vkg               =  var_Ckg(Pk_interps, lensCl_interps, nolensCl_interps, Llls, zmin, zmax, pz, bz, dnbar, fsky,\
                                     nkk=nkk,    iterative=iterative, pickleit=True, rho=rho)
      
        result            =  snr(Llls, ckg, vkg, nmodes, lmax = Lmax)
        result           /=  np.sqrt(fsky)
      
        results.append(result)

        print('\n\nFor z=%.1lf, total S/N: %.1lf to Lmax of %.1lf (fsky, thetab, DeltaT = %.3lf, %.2lf, %.2lf)' % (peakz, result, Lmax,\
                                                                                                                   fsky, thetab, DeltaT))

      results = np.array(results) 
    
      pl.semilogx(dfactors * nbar, results, ls, label=label, c=colors[ii])
  
  pl.xlabel(r'$\bar n / \rm{deg}^{2}$', fontsize=12)

  if band == 'LBG':
    pl.ylabel('Cumulative (S/N) ' + r'$/ \ \sqrt{f_{\rm{sky}}}$', fontsize=12)

  pl.xlim(1.e1, 1.e4)
  pl.ylim(0.0,  550.)
  
  pl.legend(loc=2, ncol=2)

  plt.tight_layout()
  
  pl.show()
  ##  pl.savefig('plots/%ssnr.pdf' % band)

  print('\n\nDone.\n\n')

