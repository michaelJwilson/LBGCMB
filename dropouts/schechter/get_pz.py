import  os
import  numpy as np
import  pylab as pl

from    scipy.interpolate import interp1d


def get_pz(type='g'):
    root =  os.environ['LBGCMB']
    file =  root + '/dropouts/schechter/dat/pz/pz_%s.dat' % type
    
    x, y =  np.loadtxt(file, unpack=True)

    return  interp1d(x, y, kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)


if __name__ == '__main__':
  print('\n\nWelcome.\n\n')

  dz = 0.01
  zs = np.arange(0.0, 10.0, dz)

  for type in ['BX', 'u', 'g', 'r']:
    pp   = get_pz(type=type)

    norm = np.sum(pp(zs)) * dz
        
    pl.plot(zs, pp(zs), label=type + ', %.2lf' % norm)

  pl.legend()
  pl.show()

  print('\n\nDone.\n\n')
