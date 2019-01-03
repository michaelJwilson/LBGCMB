import  numpy    as      np

from    extinct  import  extinct
from    madau    import  lephare_madau


def extinct_sed(z, ls, EBV=0.1, types=['calzetti', 'madau'], printit=False):
  """                                                                                                                                                                                                                                          Input:    rest wavelength ls [A] ordered by decreasing wavelength;    
    Returns:  Throughput array.  
  """

  throughput = np.ones_like(ls)
                                                                                                                                                                                                               
  if 'allen' in types:    
    throughput *= extinct(ls, EBV, HyperZ=True, type='allen',       atmos=False)   ## ls in angstroms                                                                                                                                      
    types.remove('allen')

  elif 'fitzpatrick' in types:
    throughput *= extinct(ls, EBV, HyperZ=True, type='fitzpatrick', atmos=False)
    types.remove('fitzpatrick')
    
  elif 'calzetti' in types:
    throughput *= extinct(ls, EBV, HyperZ=True, type='calzetti',    atmos=False)
    types.remove('calzetti')

  else:
    pass

  if 'madau' in types:
    ## exp(-tau) where tau is the optical depth due to absorption; ls in angstroms.                                                                                                                                                      
    throughput *= lephare_madau(ls, z)
    types.remove('madau')

  if len(types) != 0:
    import warnings

    print '\n\n(Warning)  Extinction rule not found for the following input: ' + ','.join(types)
        
  return throughput


if __name__ == '__main__':
  print("\n\nWelcome to extinct-sed.\n\n")

  ## Test
  z    = 2.0
  
  ls   = np.arange(1e3, 1e4, 1e2)

  tput = extinct_sed(z, ls, EBV=0.1, types=['calzetti', 'madau'], printit=False)

  print tput

  print("\n\nDone.\n\n")
  
