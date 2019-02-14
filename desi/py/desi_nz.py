import  os
import  yaml
import  numpy              as      np
import  pylab              as      pl

from    scipy.interpolate  import  interp1d   
from    desi_bz            import  lrg_bz, elg_bz, qso_bz, bgs_bz
from    utils              import  latexify


def desi_ps():
    result = {}

    for target in ['bgs_bright', 'bgs_faint', 'elg', 'lrg', 'qso']:
      if target in ['bgs_bright', 'bgs_faint']:
        fpath  = os.environ['LBGCMB'] + '/dat/desi/targets/nz_%s.dat' % 'BG'
        data   = np.loadtxt(fpath)

      else:
        fpath  = os.environ['LBGCMB'] + '/dat/desi/targets/nz_%s.dat' % target
        data   = np.loadtxt(fpath)

      dz     = data[0,1] - data[0,0]
      ngal   = data[:,2].sum()        ##  [per sq. deg].
  
      zs     = data[:,0]
      ps     = data[:,2] / ngal / dz

      result['dz_%s' % target] = dz
      result['zs_%s' % target] = zs 
      result['pz_%s' % target] = ps
   
    return  result

def get_targets():
    fpath  = os.environ['LBGCMB'] + '/dat/desi/targets/targets.yaml'
    y      = yaml.load(open(fpath))

    for target in ['elg', 'lrg', 'qso', 'bgs_bright', 'bgs_faint']:
        y['goodz_%s' % target] = y['nobs_%s' % target] * y['success_%s' % target]

    return  y 

def _get_stats():
    targets = get_targets()
    ps      = desi_ps()

    bs      = {'elg': elg_bz, 'lrg': lrg_bz, 'qso': qso_bz, 'bgs_bright': bgs_bz, 'bgs_faint': bgs_bz}

    for target in ['elg', 'lrg', 'qso', 'bgs_bright', 'bgs_faint']:
      ps['nz_%s' % target] = targets['success_%s' % target] * targets['nobs_%s' % target] * ps['pz_%s' % target]
      ps['bz_%s' % target] = bs[target] 

    return  ps

def get_stats():
    bs               =  {'elg': elg_bz, 'lrg': lrg_bz, 'qso': qso_bz, 'bgs_bright': bgs_bz, 'bgs_faint': bgs_bz}

    targets          =  get_targets()
    stats            =  _get_stats()

    ntot             =  0.0

    dz               =  0.05

    stats['dz_desi'] =    dz
    stats['zs_desi'] =  np.arange(0.0, 5.0, dz) 
    
    stats['nz_desi'] =  np.zeros_like(stats['zs_desi']) 
    stats['bz_desi'] =  np.zeros_like(stats['zs_desi'])

    for target in ['elg', 'lrg', 'qso', 'bgs_bright', 'bgs_faint']:
        ntot             +=  targets['success_%s' % target] * targets['nobs_%s' % target]  ## [per sq. deg.]

        zs                =  stats['zs_%s' % target]
        ns                =  stats['nz_%s' % target]

        ni                =  interp1d(zs, ns, kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)    
        dn                =  ni(stats['zs_desi'])
            
        ## 
        stats['nz_desi'] +=  dn
        stats['bz_desi'] +=  dn * bs[target](stats['zs_desi'])  

    ## 
    stats['bz_desi'][stats['nz_desi'] > 0.0] /= stats['nz_desi'][stats['nz_desi'] > 0.0]

    return stats
    

if __name__ == "__main__":
    import pylab as pl
    import matplotlib.pyplot as plt
    

    latexify(columns=2, equal=False, fontsize=10, ratio=None, ggplot=True, usetex=False)

    print("\n\nWelcome to the DESI p(z) calculator.\n\n")

    targets    = get_targets()
    desi_stats = get_stats()

    for target in ['elg', 'lrg', 'qso', 'bgs_bright', 'bgs_faint']:
      pl.plot(desi_stats['zs_%s' % target], np.cumsum(desi_stats['nz_%s' % target]) * desi_stats['dz_%s' % target], label=target)   

    pl.plot(desi_stats['zs_desi'], np.cumsum(desi_stats['nz_desi']) * desi_stats['dz_desi'], 'k-', label='desi')
    
    pl.xlim(3.0,  0.0)

    pl.xlabel(r'$z_{\rm{max}}$')
    pl.ylabel(r'$N(z  < z_{\rm{max}}) \ [\rm{gals. / deg}^2]$')  

    pl.legend()


    ##  Linear galaxy bias.
    ax  = pl.gca()
    ax2 = ax.twinx()

    for target in ['elg', 'lrg', 'qso', 'bgs_bright', 'bgs_faint']:
        ax2.plot(desi_stats['zs_desi'], desi_stats['bz_%s' % target](desi_stats['zs_desi']), label=target, alpha=0.6, linestyle='--')

    ax2.plot(desi_stats['zs_desi'], desi_stats['bz_desi'], 'k--', label='desi', alpha=0.6, linestyle='--')
    ax2.set_ylabel(r'$b(z)$')

    ax2.set_ylim(0.0, 4.0)

    plt.tight_layout()
    pl.show(block=True)

    print("\n\nDone.\n\n")
