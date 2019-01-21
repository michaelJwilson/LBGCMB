import numpy as np


def des_depths():
  '''
  NOTE:  SV is equivalent to the 5 yr DES depth.
  '''

  ## Assume that flux of 10 sig. source is double that of 5 sigma source:  m5 = m10 + 0.75                                                                 
  stats     = {'Y1': {'g': 24.15, 'r': 23.95, 'i': 23.25, 'z': 22.55},\
               'SV': {'g': 25.20, 'r': 25.05, 'i': 24.25, 'z': 23.65}}

  for year in np.arange(2, 6, 1):
    ## 5 sig. flux of 5 yr = sqrt(N / 5) 5 sig flux of N yr:  m_Nyr - m_SV = 2.5 log10( sqrt(N/5.) )
    stats['Y' + str(year)] = {}

    ## Add to dictionary under Y + int(year) scaling.                                                                                                
    for band in ['g', 'r', 'i', 'z']:
      stats['Y' + str(year)][band] = stats['SV'][band] + 2.5 * np.log10( np.sqrt(year / 5.) )

  return  stats


if __name__ == "__main__":
  from utils import pprint
  

  stats = des_depths()    
    
  pprint(stats)

  print("\n\nDone.\n\n")
