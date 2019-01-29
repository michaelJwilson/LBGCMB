import numpy as np
import pylab as pl

from   prep_filters import prep_filters
from   collections  import OrderedDict 
 

def get_detband(z, filters, printit=False, restframe=False):
  ##  In the case of rest frame, arbitrarily set to the r band. 
  if restframe:
    return 'r'

  ##  Observer frame. 
  le      = 1216.           ##  AA
  lo      = le * (1. + z) 

  ##  Check on the filters are ordered in wavelength.
  check   = 0.0
  
  bands   = list(filters.keys())
  inband  = np.zeros(len(bands))

  dband   = ''

  for i, band in enumerate(filters):
    mmax  = filters[band]['Ts'].max()
    hmax  = mmax / 4.

    lolim = filters[band]['ls'][filters[band]['Ts'] > hmax].min()
    hilim = filters[band]['ls'][filters[band]['Ts'] > hmax].max()

    assert (check < lolim) 

    inband[i] = (lolim <= lo) & (lo <= hilim)

    if inband[i]:
      dband = band

    if printit:
      print('%s \t %.2lf \t %.2lf \t %.2lf \t %d' % (band, lolim, hilim, lo, inband[i]))
        
    check = lolim

  ##  Check lyman-alpha is in one of the filters.
  assert np.sum(inband) > 0.0

  ##  Check it's not the last filter (in which case there 'are' non-redder (in LSST)). 
  assert dband != 'y'

  ##  Assign first detection band as the nearest redder than the filter in which lyman-alpha falls.   
  inband  = np.array([True if x == dband else False for x in bands])
  inband  = np.roll(inband, 1)

  dband   = np.array(bands)[inband]

  return  dband[0]


if __name__ == '__main__':
  print('\n\nWelcome\n\n')

  filters = prep_filters(names=['LSST', 'VIDEO'], normed=False)
  
  ##  Cut to the LSST filters.  
  olsst   = ['u', 'g', 'r', 'i', 'z', 'y']

  filters = OrderedDict(zip(olsst, [filters[x] for x in olsst]))

  zs      = np.arange(2.0, 6.5, 0.1)
  dbands  = []

  for z in zs:
    dband = get_detband(z, filters, printit=False)
    dbands.append(dband)

  mmap    = dict(zip(olsst, np.arange(6)))  
  
  dbands  = np.array([mmap[x] for x in dbands])

  for i, band in enumerate(filters.keys()):
    pl.plot(filters[band]['ls'], filters[band]['Ts'], label = filters[band]['ppkey'])

  pl.plot((1. + zs) * 1216., dbands)
  pl.legend()
  pl.show(block=True)
  
  print('\n\nDone.\n\n')
