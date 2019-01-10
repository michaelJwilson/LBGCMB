import  os
import  numpy            as      np
import  pylab            as      pl

from    bz               import  get_dropoutbz
from    prep_camb        import  CAMB
from    pmh              import  Pmm, get_PkInterps
from    params           import  get_params
from    cosmo            import  cosmo
from    scipy.integrate  import  nquad
from    growth_rate      import  growth_rate 
from    utils            import  latexify


if __name__ == '__main__':
    print('\n\nWelcome to the RSD S/N calculator.')

    latexify(columns=1, equal=True, fontsize=10, ratio=None, ggplot=True, usetex=True)

    ##  drop_0.2000_12.60_0.60.pxs  
    ##  drop_0.2000_12.83_0.45.pxs  
    ##  drop_0.2500_12.33_0.40.pxs

    root = os.environ['LBGCMB']
    dat  = np.loadtxt(root + '/dropouts/Martin/dat/summary_v2/drop_0.2500_12.33_0.40.pxs')

    pl.loglog(dat[:,0],        dat[:,1],  'b-', label=r'$z=3$', alpha=0.8)
    pl.loglog(dat[:,0], np.abs(dat[:,2]), 'b-', label='', alpha=0.8)

    root = os.environ['LBGCMB']
    dat  = np.loadtxt(root + '/dropouts/Martin/dat/summary_v2/drop_0.2000_12.60_0.60.pxs')

    pl.loglog(dat[:,0],        dat[:,1],  'g-', label=r'$z=4$', alpha=0.5)
    pl.loglog(dat[:,0], np.abs(dat[:,2]), 'g-', label='', alpha=0.5)
    
    pl.xlabel(r'$k_{\rm{max}} \ [(h \ \rm{Mpc}^{-1})]$')
    pl.ylabel(r'$P_{\ell}(k) \ [(h^{-1} \rm{Mpc})^3]$')

    pl.xlim(2.e-2, 0.5)
    pl.ylim(1.e2, 2.e4)

    pl.legend(loc=1)

    pl.savefig('plots/rsd_multipoles.pdf', bbox_inches='tight')

    print('\n\nDone.\n\n')
