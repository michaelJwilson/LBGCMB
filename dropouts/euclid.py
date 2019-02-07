import numpy as np
import pylab as pl

from scipy.interpolate import interp1d

def euclid_bz(z):
    return  np.sqrt(1. + z)

def euclid_area():
    return  15000.          ##  [deg2]

@np.vectorize
def euclid_nz(z):
    zs = np.array([1.4, 1.5, 1.8, 2.0])
    ns = np.array([10., 8.0, 3.5, 1.1]) * 1.e-4

    it = interp1d(zs, ns, kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)

    return  it(z)


if __name__ == '__main__':
    zs = np.arange(0.0, 2.5, 0.01)
    ns = euclid_nz(zs)

    pl.plot(zs, ns)
    pl.show()
    
    print('\n\nDone.\n\n')
