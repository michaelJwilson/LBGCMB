import  os 
import  numpy             as      np
import  pylab             as      pl

from    reddy.specs       import  samplestats
from    get_schechters    import  get_schechters
from    nbar              import  comovdensity
from    dVols             import  dVols
from    cosmo             import  cosmo
from    params            import  get_params
from    nbar              import  dndz
from    gen_pz            import  get_pz
from    reddy.pz          import  get_pz as reddy_getpz
from    scipy.signal      import  medfilt
from    scipy.interpolate import  interp1d

params = get_params()

def get_completeness(plotit=False, save=False, interp=True):
    mlim           =  25.5
    key            = 'BX'
    stats          =  samplestats(printit = False)

    midz, alpha, M_star, phi_star = get_schechters(stats, key)
    
    zs             =  np.arange(0.0, 10.0, 0.01)
    zs, dVs, ns    =  dndz(zs, phi_star, M_star, alpha, mlim, type='app', printit = True, completeness=None, app_linelim=False)
    zs, ps         =  get_pz(zs, ns, lambda z:  1.0)

    red_pz         =  reddy_getpz(interp=True)

    Cs             =  red_pz(zs) / ps
    Cs            /=  Cs.max()
    Cs             =  medfilt(Cs, kernel_size=11)
    
    if plotit:
      pl.plot(zs, Cs)

      pl.legend(frameon=False)
      pl.show()

    if save:
      ##  And save. 
      np.savetxt(os.environ['LBGCMB'] + 'dat/completeness.dat')

    if interp:
        return  interp1d(zs, Cs, kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)

    else:
        return  zs, Cs

if __name__ == '__main__':
    zs, ps  = get_completeness(plotit=True, save=False, interp=False)
    
    print('\n\nDone.\n\n')
