import  os
import  yaml
import  numpy              as      np
import  pylab              as      pl

from    scipy.interpolate  import  interp1d   
from    desi_bz            import  lrg_bz, elg_bz, qso_bz, bgs_bz


def desi_ps():
    result = {}

    for target in ['BG', 'elg', 'lrg', 'qso']:
      fpath  = os.environ['LBGCMB'] + '/dat/desi/targets/nz_%s.dat' % target
      data   = np.loadtxt(fpath)

      dz     = data[0,1] - data[0,0]
      ngal   = data[:,2].sum()        ##  [per sq. deg].
  
      zs     = data[:,0]
      ps     = data[:,2] / ngal / dz

      result['zs']             = zs 
      result['pz_%s' % target] = ps
   
    return  result

def get_targets():
    fpath  = os.environ['LBGCMB'] + '/dat/desi/targets/targets.yaml'
    y      = yaml.load(open(fpath))

    for target in ['elg', 'lrg', 'qso', 'bgs_bright', 'bgs_faint']:
        y['goodz_%s' % target] = y['nobs_%s' % target] * y['success_%s' % target]

    return  y 

def get_stats():
    targets = get_targets()
    ps      = desi_ps()

    bs      = {'elg': elg_bz, 'lrg': lrg_bz, 'qso': qso_bz, 'bgs_bright': bgs_bz, 'bgs_faint': bgs_bz}

    for target in ['elg', 'lrg', 'qso', 'bgs_bright', 'bgs_faint']:
      if target in ['bgs_bright', 'bgs_faint']:
          ps['nz_%s' % target] = targets['goodz_%s' % target] * ps['pz_%s' % 'BG']

      else:
          ps['nz_%s' % target] = targets['goodz_%s' % target] * ps['pz_%s' % target]

      ps['bz_%s' % target]     = bs[target] 

    print(ps)


if __name__ == "__main__":
    print("\n\nWelcome to the DESI p(z) calculator.\n\n")

    desi_stats = get_stats()

    print("\n\nDone.\n\n")
