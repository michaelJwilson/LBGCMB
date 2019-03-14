import  os
import  numpy              as      np

from    scipy.interpolate  import  interp1d


def get_completeness():
  root           =  os.environ['LBGCMB']
  path           =  root + '/dropouts/Malkan/dat/completeness.dat'

  completeness   =  np.loadtxt(path)
  interp         =  interp1d(completeness[:,0], completeness[:,1], kind='linear', assume_sorted=False, bounds_error=False, fill_value=0.0)

  return  interp
