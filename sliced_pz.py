import  numpy  as  np


def sliced_pz(zz, zmin, zmax, survey_pz):
  '''
  Return p(zz) for a normalisation forced to be unity over 
  [zmin, zmax], rather than 0 < z < 10.
  '''

  zs                = np.arange(zmin, zmax, 0.01)
  norm              =  np.trapz(survey_pz(zs), zs)

  result            = survey_pz(zz) / norm

  result[zz < zmin] = 0.0
  result[zz > zmax] = 0.0

  return  result

