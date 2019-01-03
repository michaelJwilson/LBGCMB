import  os
import  numpy              as      np
import  pylab              as      pl

from    scipy.interpolate  import  interp1d   


root  = os.environ['LBGCMB']
root += '/desi/dat/'

def elg_pz(z, type=''):
    ##  ['nz_elg.dat', 'nz_elg_cdr.dat', 'nz_elg_fdr.dat']
    fpath  = root + 'nz_elg%s.dat' % type
    data   = np.loadtxt(fpath)

    dz     = data[1,0] - data[0,0]
    norm   = data[:,2].sum() 
  
    zee    = data[:,0]
    pz     = data[:,2] / norm / dz

    interp = interp1d(zee, pz, kind='cubic', bounds_error=False, fill_value=0.0)

    return interp(z)

if __name__ == "__main__":
    print("\n\nWelcome to the ELG p(z) calculator.")

    zs = np.arange(0.0, 4.0, 0.1)

    elg_pz(zs, type='')

    print("\n\nDone.\n\n")
