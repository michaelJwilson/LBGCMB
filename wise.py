import  os
import  numpy              as      np
import  pylab              as      pl

from    scipy.integrate    import  simps
from    scipy.interpolate  import  interp1d


def get_wisesn():
    return  1.2e-10

def get_wisepz():
    root = os.environ['LBGCMB']

    data = np.loadtxt(root + '/dat/wise.dat')
    
    z    = data[:,0]
    dz   = data[1,0] - data[0,0]

    pz   = data[:,1]

    norm = dz * np.sum(pz)
    pz  /= norm

    pz   = interp1d(z, pz, kind='linear', copy=True, fill_value=0.0, assume_sorted=False, bounds_error=False)

    return  pz
    

if __name__ == '__main__':
    print('\n\nWelcome to WISE p(z).\n\n')
    
    pz   = get_wisepz()

    zs   = np.arange(0.0, 4.0, 0.01)

    pl.plot(zs, pz(zs))

    pl.show()
    
    print('\n\nDone.\n\n')
