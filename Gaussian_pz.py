import  numpy              as      np
from    scipy.integrate    import  simps


def Gaussian(z, z0, sigma):
  arg  = -(z - z0)**2.
  arg /=   2.*sigma**2.

  return np.exp(arg)

def Gaussian_pz(z, z0 = 2.96, sigma = 0.24):
  """
  Gaussian form for p(z) centered on z0 with standard
  deviation of sigma;  Note:  integral is normalised 
  to unity.

  Note:  defaults to Hildbrandt (2009) u-drops. 
  """

  ## Normalise \int p(z) dz to be unity for zmin = 0 and zmax = 10.                                                                                        
  zs      = np.arange(0., 10., 0.01)

  ## Integrate Gaussian p(z) to get the normalisation.                                                                                                     
  norm    = simps(Gaussian(zs, z0, sigma), dx=0.01)

  ## p(z) normalised to unity over 0 < z < 10.                                                                                                             
  return  Gaussian(z, z0, sigma)/norm


if __name__ == "__main__":
  print("\n\nWelcome to Gaussian p(z) generator.\n\n")

  z0   = 3.0 
  nbar = 100.  ## [deg2] 


  zs   = np.arange(1.0, 4.0, 0.1)

  for zee in zs:
    print  Gaussian_pz(zee, z0, nbar)

  print("\n\nDone.\n\n")
