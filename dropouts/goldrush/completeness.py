import  os
import  numpy              as      np

from    scipy.interpolate  import  interp1d 


def get_completeness(drop='g'):
  mmap           =  {'g': 4, 'r': 5}

  root           =  os.environ['LBGCMB']
  path           =  root + '/dropouts/goldrush/cats/completeness/completeness_z%d.dat' % mmap[drop]

  print('\nReading:  %s' % path)

  completeness   =  np.loadtxt(path)
  interp         =  interp1d(completeness[:,0], completeness[:,1], kind='linear', assume_sorted=False, bounds_error=False, fill_value=0.0)

  return  interp


if __name__ == "__main__":
  print("\n\nWelcome to completeness.\n\n")

  get_completeness('r')

  print("\n\nDone.\n\n")
