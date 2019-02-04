import  numpy              as      np
import  pylab              as      pl

from    prep_camb          import  CAMB
from    prep_Llls          import  prep_Llls
from    pmh                import  Pmm, get_PkInterps, linz_bz
from    Ckg                import  Ckg 
from    Cgg                import  Cgg, Ngg
from    Gaussian_pz        import  Gaussian_pz
from    bz                 import  get_dropoutbz
from    completeness       import  get_dropoutpz
from    scipy.interpolate  import  interp1d
from    cib                import  nWCIB
from    lensing            import  Ckk, var_Ckk
from    wise               import  get_wisepz
from    Nkk                import  Nkk
from    bolometers         import  bolometers
from    utils              import  latexify
from    whitebook_pz       import  const_pz 
from    schmittfull_nz     import  get_ss17_samples


latexify(fig_width=None, fig_height=None, columns=1, equal=True, fontsize=10)

cmbexp                             = 'CMBS4'
plotit                             =   True

cambx                              =  CAMB()
Pk_interps                         =  get_PkInterps(cambx)

(lensCl_interps, nolensCl_interps) =  cambx.get_Cls()

NLlls, Llls, nmodes                =  prep_Llls(NLlls = 60, Lmin = 50., Lmax = 5000., log10=True)

fsky, thetab, DeltaT, iterative    =  bolometers[cmbexp]['fsky'],   bolometers[cmbexp]['thetab'],\
                                      bolometers[cmbexp]['DeltaT'], bolometers[cmbexp]['iterative']

ns, ps, bs, ss                     =  get_ss17_samples(nolsst=False)

samples                            =  zip(bs, ps, ns)
samples                            =  samples[2:4]

print(samples)

exit(1)

result                             =  []

## Ckk
ckk                                =  Ckk(Pk_interps, Llls, pickle=False)
nkk                                =  Nkk(lensCl_interps, nolensCl_interps, Llls, terms=['TT', 'TE', 'EE', 'EB'],\
                                          thetab=thetab, DeltaT=DeltaT, iterative=iterative, pickleit=True)

Llls                               =  Llls[::5]
ckk                                =   ckk[::5]
nkk                                =   nkk[::5]

zmin                               =   0.01
zmax                               =  10.00

for ii, L in enumerate(Llls):
  L              =  np.array([L])
  kk             =  ckk[ii]
  
  kg             = [ Ckg(Pk_interps, L, zmin, zmax, p, b, zeff = False)                           for [b, p, n] in samples]
  
  ## Shotnoise to be added?
  gg             = [[Cgg(Pk_interps, L, zmin, zmax, p, b, bz2  = b2, survey_pz2 = p2, zeff=False) for [b, p, n] in samples] for [b2, p2, n2] in samples]
  
  kg             =  np.array(kg)
  gg             =  np.array(gg)[:,:,0]

  for i, [b, p, n] in enumerate(samples):
    gg[i,i]     +=  Ngg(L, zmin, zmax, p, n)

  diag           =  np.diag(gg)
  rho            =  np.copy(gg)

  for i in np.arange(len(diag)):
    for j in np.arange(len(diag)):
      rho[i,j]  /=  np.sqrt(diag[i]) * np.sqrt(diag[j])
  
  kg            /=  np.sqrt(kk) * np.sqrt(diag)[:, None]

  irho           =  np.linalg.inv(rho)
  interim        =  0.0

  for i in np.arange(len(kg)):
    for j in np.arange(len(kg)):
      interim   += kg[i] * irho[i, j] * kg[j] 
      
  result.append([L, np.sqrt(interim)])
  

result  = np.array(result)[:,:,0]
rho     = interp1d(result[:,0], result[:,1], kind='linear', bounds_error=False, fill_value=0.0, assume_sorted=False)

np.savetxt('rho/withlsst.txt', np.c_[Llls, (1. - rho(Llls) ** 2.)])

## pl.loglog(Llls, Llls * (1. - rho(Llls) ** 2.) * ckk, label=r'$L \cdot (1 - \rho_L^2) \ C_{\kappa \kappa}$')

pl.loglog(Llls, Llls * ckk,      label=r'$L \cdot C_{\kappa \kappa}$')
pl.loglog(Llls, Llls * nkk, 'k', label=r'$L \cdot N_{\kappa \kappa}$', alpha=0.5, dashes=[3,1])

if plotit:
  for type in ['no-lsst', 'lsst']:
    data  = np.loadtxt('rho/%s.txt' % type)
    label = '%s ' % type + '$L \cdot C_{\kappa \kappa}$'

    pl.loglog(data[:,0], data[:,0] * data[:,1] * ckk, label=label)

pl.xlabel(r'$L$')

pl.xlim(50.,    4.e3)
pl.ylim(1.e-6, 3.e-5)

pl.legend(ncol=2, handlelength=.5, loc=3, handletextpad=0.5)

pl.savefig('plots/yu.pdf', bbox_inches='tight')

print('\n\nDone.\n\n')
