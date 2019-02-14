import  os
import  yaml
import  numpy              as      np
import  pylab              as      pl

from    scipy.interpolate  import  interp1d   


def pz(z, type=''):
    

    ##  [nz_BG.dat, nz_elg.dat, nz_elg_cdr.dat, nz_elg_fdr.dat, nz_lrg.dat, nz_lrg_cdr.dat, nz_qso.dat]
    fpath  = os.environ['LBGCMB'] + '/dat/desi/targets/nz_elg%s.dat' % type
    data   = np.loadtxt(fpath)

    dz     = data[1,0] - data[0,0]
    norm   = data[:,2].sum() 
  
    zee    = data[:,0]
    pz     = data[:,2] / norm / dz

    interp = interp1d(zee, pz, kind='cubic', bounds_error=False, fill_value=0.0)

    return interp(z)

def get_targets():
    fpath  = os.environ['LBGCMB'] + '/dat/desi/targets/targets.yaml'
    y      = yaml.load(open(fpath))

    for target in ['elg', 'lrg', 'qso', 'bgs_bright', 'bgs_faint']:
        y['goodz_%s' % target] = y['nobs_%s' % target] * y['success_%s' % target]

        print(target, y['goodz_%s' % target])

if __name__ == "__main__":
    print("\n\nWelcome to the DESI p(z) calculator.")

    ##  zs = np.arange(0.0, 4.0, 0.1)
    ##  ps = pz(zs, type='')

    get_targets()

    print("\n\nDone.\n\n")
