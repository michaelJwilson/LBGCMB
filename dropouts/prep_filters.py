import  collections

import  os
import  numpy             as np
import  astropy.constants as const


def prep_filters(filter='LSST'):
  filters  = collections.OrderedDict()                                              ## Generate the required filters; Euclid: Y, J and H.                    

  if filter == 'LSST':
    for band in ['u', 'g', 'r', 'i', 'z', 'y']:
      ppkey  = r"LSST-$" + "%s" % band + r"$"                                       ## Pretty print version of key.  

      fname  = os.environ['LBGCMB'] + '/dropouts/filters/lsst/total_%s.dat' % band

      ## Check transmission in energy or number;  lsst filters in nanometres.
      ## Lsst throughput (photon counting) in wavelength; [0 to 1].
      ## LSST:  nanometers on input to Angstroms.
      data   = np.loadtxt(fname)

      ## Filters defined in wavelength; conversion for integral over frequency.
      ls     = 10. * data[:,0]
      vs     = (1e10/ls) * const.c.to('m/s').value

      filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1]}
  
  elif filter == 'VIDEO':
    for band in ['Y', 'J', 'H', 'K']:
      ppkey  = r"VIDEO-$" + "%s" % band + r"$"                                       ## Pretty print version of key.

      fname  = "filters/video/%s.txt" % band
      data   = np.loadtxt(fname)

      ls     = data[:,0]
      vs     = (1e10/ls) * const.c.to('m/s').value

      filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1] / data[:,1].max()}

  else:
    raise ValueError('Requested filter set is not available.')
  
  return  filters
  

if __name__ == "__main__":
  import pylab as pl

  print("\n\nWelcome to prep-filters.\n\n")

  filters = prep_filters()

  for i, band in enumerate(filters.keys()):
    print band

    pl.plot(filters[band]['ls'], filters[band]['Ts'], label = filters[band]['ppkey'])

  ## for filter in filters.values():
  ##  print filter['ppkey']

  pl.legend(ncol=2)

  pl.savefig('plots/lsst-filters.pdf')

  print("\n\nDone.\n\n")
