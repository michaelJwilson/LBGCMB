import  numpy  as  np


def  sliced_pz(zz, zmin, zmax, survey_pz):
  '''
  Return p(zz) for a normalisation forced to be unity over 
  [zmin, zmax], rather than 0 < z < 10.
  '''

  dz                = 0.001
  zs                = np.arange(zmin, zmax, dz)
  ps                = survey_pz(zs)

  ##  Normalisation between zmin and zmax.
  norm              = np.trapz(ps, zs)
  result            = survey_pz(zz) / norm

  result[zz < zmin] = 0.0
  result[zz > zmax] = 0.0

  ##  Constant linear spacing of z array (to six deciaml places). 
  assert  len(set(np.around(np.ediff1d(zz), 6))) == 1
  
  dz      = list(set(np.around(np.ediff1d(zz), 6)))[0]
  result /= (np.sum(result) * dz)

  return  result


if __name__ == '__main__':
  from  Gaussian_pz  import  Gaussian_pz


  print('\n\nWelcome.\n\n')

  zmin   = 3. 
  zmax   = 4.

  dz     = 0.05
  zz     = np.arange(3., 4., dz)
  pz     = Gaussian_pz

  result = sliced_pz(zz, zmin, zmax, pz)

  print(np.sum(result) * dz)
  
  print('\n\nDone.\n\n')
  
