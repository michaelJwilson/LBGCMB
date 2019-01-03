import  numpy              as      np
import  pylab              as      pl

from    prep_camb          import  CAMB
from    prep_Llls          import  prep_Llls
from    pmh                import  Pmm, get_PkInterps, linz_bz
from    Ckg                import  Ckg 
from    Cgg                import  Cgg
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


latexify(fig_width=None, fig_height=None, columns=1, equal=True, fontsize=10)

cmbexp                             = 'CMBS4'

cambx                              =  CAMB()
Pk_interps                         =  get_PkInterps(cambx)

(lensCl_interps, nolensCl_interps) =  cambx.get_Cls()

zee, pzee                          =  get_dropoutpz()
pz                                 =  interp1d(zee, pzee, kind='linear', bounds_error=False, fill_value=0.0, assume_sorted=False)

fsky, thetab, DeltaT, iterative    =  bolometers[cmbexp]['fsky'],   bolometers[cmbexp]['thetab'],\
                                      bolometers[cmbexp]['DeltaT'], bolometers[cmbexp]['iterative']

## dropout_bz                      =  get_dropoutbz()
wise_pz                            =  get_wisepz()

bs                                 =  [linz_bz,     linz_bz]
ps                                 =  [pz,      Gaussian_pz]

samples                            =  zip(bs, ps)
result                             =  []

NLlls, Llls, nmodes                =  prep_Llls(NLlls = 60, Lmin = 50., Lmax = 5000., log10=True)

## Ckk
ckk                                =  Ckk(Pk_interps, Llls, pickle=False)
nkk                                =  Nkk(lensCl_interps, nolensCl_interps, Llls, terms=['TT', 'TE', 'EE', 'EB'],\
                                          thetab=thetab, DeltaT=DeltaT, iterative=iterative, pickleit=True)

ckk                                = ckk[::5]
nkk                                = nkk[::5]

pl.loglog(Llls[::5], Llls[::5] * ckk,      label=r'$L \cdot C_{\kappa \kappa}$')
pl.loglog(Llls[::5], Llls[::5] * nkk, 'k', label=r'$L \cdot N_{\kappa \kappa}$', alpha=0.5, dashes=[3,1])

for zmax in 0.1 * np.arange(2., 31., 2.):
  pl.loglog(Llls, Llls * Ckk(Pk_interps, Llls, pickle=False, zmax=zmax), c='k', lw=0.1)

for zmax in [2., 3., 4., 5., 7., 9.]:
  pl.loglog(Llls, Llls * Ckk(Pk_interps, Llls, pickle=False, zmax=zmax), label=r'$z=%.0lf$' % zmax, lw=1.0)

## pl.loglog(Llls[::5], (1. - rho(Llls[::5]) ** 2.) * ckk, label=r'$(1 - \rho_L^2) \ C_{\kappa \kappa}$')

pl.xlabel(r'$L$')

pl.xlim(50.,    4.e3)
pl.ylim(1.e-7, 2.e-4)

pl.legend(ncol=3, handlelength=.5, loc=1, handletextpad=0.5)

pl.savefig('plots/ckkz.pdf', bbox_inches='tight')

