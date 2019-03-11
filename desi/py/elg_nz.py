import os 
import numpy as np
import pylab as pl

from   scipy.interpolate import interp1d


def get_pz(dz=0.1, interp=False):
    root  = os.environ['LBGCMB'] + '/desi/dat/'
    fpath = root + '/nz_elg.dat'

    dat   = np.loadtxt(fpath)
    ngal  = np.sum(dat[:,2]) 
    
    zs    = 0.5 * (dat[:,0] + dat[:,1])
    pz    = dat[:,2] / ngal / dz
    
    if interp:
        return  interp1d(zs, pz, kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)

    else:
        return  pz


if __name__ == '__main__':
    print('\n\nWelcome to ELG p(z).\n\n')

    zs  = np.arange(0.0, 2.0, 0.1)
    ps  = get_pz()

    pl.plot(zs, ps)
    pl.show()

    print('\n\nDone.\n\n')
