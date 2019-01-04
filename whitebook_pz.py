import numpy as np


def const_pz(z, zlo=0.0, zhi=1.0):
    if isinstance(z, float):
        if (zlo < z) & (z < zhi):
            return  1. / (zhi - zlo)

        else:
            return  0.0
    else:
        z      = np.asarray(z)
        result = np.ones_like(z) / (zhi - zlo)

        result[z < zlo] = 0.0
        result[z > zhi] = 0.0
    
        return  result

def whitebook_pz(z, ilim = 25.3):
  z0       = 0.0417*ilim - 0.744                    ##  LSST science book (3.8);                                                                                      

  ## e.g.  46 galaxies per arcmin^2 for ilim = 25.
  ## Ng    = 46.*10.**(0.31*(ilim - 25.))           ##  LSST science book (3.7);                                                                                      

  zz       = z / z0
  result   = 0.5 * zz ** 2. * np.exp(-zz) / z0      ##  integral dz from 0 to inf is unity.                                                                            

  return  result

def whitebook_dpdz(z, ilim):
  z0       = 0.0417*ilim - 0.744                    ##  LSST science book (3.8); 

  zz       = z/z0
  result   = 0.5*(2. * zz - zz ** 2.) * np.exp(-zz) / (z0 ** 2.)

  return  result

def Ng(ilim, fmask=0.12):
  '''
  ## LSST science book (deprecated; e.g. Modi ++).
  result  = 46.0 * 10. ** (0.31 * (ilim - 25.0))
  result *= 60. * 60.
  '''

  ## Latest: https://arxiv.org/pdf/1809.01669.pdf
  ## We define Y1:  24.1 and Y10:  25.3
  result  = 42.9 * (1. - fmask) * 10. ** (0.359 * (ilim - 25.0))
  result *= 60. * 60.

  return  result

def lsst_chang_pz(z, ilim=24.1, source=False):
  dz = 0.01
  zs = np.arange(0., 10., dz)

  if source is False:
    ##  Distributions for clustering sources. 
    if ilim == 24.1:
        ##  Y1
        lengthscale = 0.26
        exponent    = 0.94

    elif ilim == 25.3:
        ##  Y10
        lengthscale = 0.28
        exponent    = 0.90

    else:
        raise ValueError("\n\nRequested ilim not available for LSST (Chang) p(z). \n\n")

  elif source is True:
    ##  Distribution for shear sources.
    if ilim == 24.1:
        ##  Y1                                                                                                
        lengthscale = 0.13
        exponent    = 0.78

    elif ilim == 25.3:
        ##  Y10                                                                                              
        lengthscale = 0.11
        exponent    = 0.68

  norm    = zs * zs * np.exp(- (zs / lengthscale) ** exponent )
  result  =  z *  z * np.exp(- (z  / lengthscale) ** exponent )

  result /= np.sum(norm * dz)

  return  result


if __name__ == '__main__':
  import pylab as pl
  import numpy as np


  print('\n\nWelcome to whitebook pz.\n\n')
  
  dz       = 0.01
  zs       = np.arange(0.0, 2.51, dz)
  '''
  for ilim in [21., 23., 25.3]:
    Ng     = 1.7 * 10. ** (5. + 0.31 * (ilim - 25.)) ## [deg^-2], McQuinn and White eqn. (3).

    pl.semilogy(zs, Ng * whitebook_pz(zs, ilim))
    pl.semilogy(zs, Ng * whitebook_dpdz(zs, ilim), '--')

  pl.ylim(1., 3.e5)

  pl.show()
  '''

  print const_pz(zs, zlo=2.0, zhi=2.3)

  print('\n\nDone.\n\n')
