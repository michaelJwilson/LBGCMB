import  numpy   as      np

from    params  import  get_params


params = get_params()

def DetectorNoise(ell, thetab = 1., DeltaT = 1., type='TT'):
  """
  Return the expected instrumental noise of the CMB experiment. 
  """

  DeltaT        *= 10.**-6.              ##  Convert from micro K to K.                                                            

  arcmin2radians = np.pi/(60.*180.)

  DeltaT        *= arcmin2radians        ##  Convert from K-arcmin to K-radian.                                                    
  thetab        *= arcmin2radians        ##  Arcmin to radians.                                                                 

  DeltaP         = np.sqrt(2.) * DeltaT  ##  Fully polarised detectors.                                                            
  
  if type == 'TT':
    return (DeltaT/params['Tcmb'])**2.*np.exp(ell*(ell + 1.)*thetab**2./(8.*np.log(2.))) ##  DeltaT in K/radian, converted from muK/arcmin.                  

  elif type in ['EE', 'BB', 'EB']:
    return (DeltaP/params['Tcmb'])**2.*np.exp(ell*(ell + 1.)*thetab**2./(8.*np.log(2.))) ##  DeltaP = root(2) * DeltaT                                      

  elif type == 'TE':
    return 0.0                           ## No instument noise added to TE                                                                                   

  else:
    raise ValueError("Noise estimate is not available for %s combination." % type)


if __name__ == '__main__':
  from  prep_Llls import prep_Llls


  print('\n\nWelcome to Detector noise.\n\n')

  Llls, nmodes = prep_Llls()

  print  DetectorNoise(Llls, thetab = 1., DeltaT = 1., type='TT')

  print('\n\nDone.\n\n')
