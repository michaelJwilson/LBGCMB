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
from    pmh                      import  Pmm, get_PkInterps, linz_bz


print("\n\nWelcome to prep. Cgg.\n\n")

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

print('\n\nDone.\n\n')
