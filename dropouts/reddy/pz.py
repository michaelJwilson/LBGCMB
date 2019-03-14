import os
import numpy as np

from   scipy.interpolate import interp1d


def get_nz():
  zs, ns    =  np.loadtxt(os.environ['LBGCMB'] + '/dropouts/reddy/dat/nz.dat', unpack=True)

  return  zs, ns

def get_pz(interp = True):
  zs, ns    =   get_nz()

  ngal      =  ns.sum()

  dz        =  zs[1] - zs[0]

  ##  Fraction of galaxies in that bin.                                                                                                                  
  ps        =  ns / ngal

  ##  Probability density.                                                                                                                               
  ps       /=  dz

  if interp:
    return  interp1d(zs, ps, kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)

  else:
    return  zs, ps


if __name__ == '__main__':
    zs, ps  = get_pz(interp = False) 
