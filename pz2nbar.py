import  numpy        as      np

from    Gaussian_pz  import  Gaussian_pz


def pz_slice(zmin, zmax, pz):
    '''                                                                                                                                                    
    Returns pz integrated between zmin and zmax.
    Relies on p(z) being correctly normalised. 
    '''
    dz      = 0.001
    zz      = np.arange(zmin, zmax, dz)

    ##  Integrate over unnormalised survey_pz from zmin to zmax.                                                                                            
    result  = np.trapz(ps(zz), zz, dx=dz)

    return  result

def nbar_convert(nbar, unit='arcmin'):                                                                                                                     
  if unit == 'arcmin':                                                                                                                                  
      return nbar/60./60.                                                                                                                         
  
  elif unit == 'str':                                                                                                                                      
      ## One degree is np.pi/180. radians                                                                                                                 
      deg_str = np.pi/180.                                                                                                                                  

      ## Galaxies per steradian.                                                                                                                             
      return nbar/deg_str**2.                                                                                                                        
  
  else:                                                                                                                                                   
    raise ValueError("Requested unit is not available for dropout p(z).")                                                                                  


if __name__ == "__main__":
  from  Gaussian_pz  import  Gaussian_pz

  print("\n\nWelcome to Gaussian p(z) generator.\n\n")

  zmin     = 2.0
  zmax     = 3.0                                                                                                                                   
  nbar     = 27.0  ## [deg2]


  result   = pz_slice(zmin, zmax, Gaussian_pz)
  result  *= nbar

  print  nbar_convert(nbar, unit='str')

  print("\n\nDone.\n\n")    
