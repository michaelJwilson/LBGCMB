import os 
import numpy as np
import pylab as pl

def get_wschechter(z):
  root                        = os.environ['LBGCMB']
  fpath                       = root + '/dropouts/Martin/dat/schechter/schechters.dat'

  zs, MUVs, phi_stars, alphas = np.loadtxt(fpath, unpack=True)
  phi_stars                  *= 1.e-3

  index                       = np.abs(zs - z) == np.abs(zs - z).min()
  midz                        = zs[index][0] 
  alpha                       = alphas[index][0]
  MUV                         = MUVs[index][0]
  phi_star                    = phi_stars[index][0]
  
  return  midz, alpha, MUV, phi_star


if __name__ == '__main__':
    print('\n\nWelcome to White Schechters.\n\n')

    midz, alpha, MUV, phi_star = get_wschechter(z = 4.7)

    print(midz)
    print(alpha)
    print(MUV)
    print(phi_star)

    print('\n\nDone.\n\n')
