import  os 
import  collections

import  numpy             as np
import  astropy.constants as const


def prep_filters(names=['LSST', 'VIDEO'], normed=False):
  filters  =  collections.OrderedDict()                                              ## Generate the required filters; Euclid: Y, J and H.                  
  root     =  os.environ['THREEDHST']

  if 'LSST' in names:
    print('\n\nAdding LSST filters (u, g, r, i, z and y).')

    for band in ['u', 'g', 'r', 'i', 'z', 'y']:
      ppkey  = r"LSST-$" + "%s" % band + r"$"                                       ## Pretty print version of key.  

      fname  = root + "/filters/lsst/total_%s.dat" % band

      ## Check transmission in energy or number;  lsst filters in nanometres.
      ## Lsst throughput (photon counting) in wavelength; [0 to 1].
      ## LSST:  nanometers on input to Angstroms.
      data   = np.loadtxt(fname)

      ## Filters defined in wavelength; conversion for integral over frequency.
      ls     = 10. * data[:,0]
      vs     = (1e10/ls) * const.c.to('m/s').value

      if normed:
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1] / data[:,1].max()}
      
      else:
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1]}
        
    names.remove('LSST')

  if 'VIDEO' in names:
    print('Adding VIDEO filters (Y, J, H and K).')

    for band in ['Y', 'J', 'H', 'K']:
      ppkey  = r"VIDEO-$" + "%s" % band + r"$"                                       ## Pretty print version of key.

      fname  = root + "/filters/video/%s.txt" % band
      data   = np.loadtxt(fname)

      ls     = data[:,0]
      vs     = (1e10/ls) * const.c.to('m/s').value

      if normed:
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1] / data[:,1].max()}

      else:   
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1]}

    names.remove('VIDEO')

  if 'STEIDEL' in names:
    print('Adding Steidel filters (U, G and R).')

    for band in ['U', 'G', 'R']:
      ppkey  = r"STEIDEL-$" + "%s" % band + r"$"
      fname  = root + "/filters/steidel/%s.pb" % band
      data   = np.loadtxt(fname)

      ls     = data[:,0]
      vs     = (1.e10 / ls) * const.c.to('m/s').value

      if normed:
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1] / data[:,1].max()}

      else:   
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1]}

    names.remove('STEIDEL')

  if 'SUBARU' in names:
    print('Adding Subaru filters (B and V).')

    for band in ['B', 'V']:
      ppkey  = r"SUBARU-$" + "%s" % band + r"$"
      fname  = root + "/filters/subaru/%s.pb" % band
      data   = np.loadtxt(fname)

      ls     = data[:,0]
      vs     = (1.e10 / ls) * const.c.to('m/s').value

      if normed:
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1] / data[:,1].max()}

      else:   
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1]}

    names.remove('SUBARU')

  if 'JKC' in names:
    print('Adding JKC filters (I).')

    for band in ['I']:
      ppkey  = r"JKC-$" + "%s" % band + r"$"
      fname  = root + "/filters/jkc/%s.pb" % band
      data   = np.loadtxt(fname)

      ls     = data[:,0]
      vs     = (1.e10 / ls) * const.c.to('m/s').value

      if normed:
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1] / data[:,1].max()}

      else:   
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1]}

    names.remove('JKC')

  if 'HUBBLE' in names:
    print('Adding HUBBLE filters (acs_f435w, acs_f606w, acs_f775w, acs_f850lp).')

    for band in ['acs_f435w', 'acs_f606w', 'acs_f775w', 'acs_f850lp']:
      ppkey  = r"HUBBLE-$" + "%s" % band.upper() + r"$"
      fname  = root + "/filters/hst/%s.pb" % band
      data   = np.loadtxt(fname)

      ls     = data[:,0]
      vs     = (1.e10 / ls) * const.c.to('m/s').value

      print(band.upper())

      if normed:
        filters[band.upper()] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1] / data[:,1].max()}

      else:
        filters[band.upper()] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1]}

    names.remove('HUBBLE')

  if len(names) != 0:
    ##  Catch any requested filters that are not either LSST or VIDEO.
    raise ValueWarning('Requested filters are not available:' + '  '.join(x for x in name))

  return  filters
  

if __name__ == "__main__":
  import pylab as pl


  print("\n\nWelcome to prep-filters.\n\n")

  ##  filters = prep_filters(['LSST', 'VIDEO', 'STEIDEL', 'SUBARU', 'JKC', 'HUBBLE'], normed=True)
  filters = prep_filters(['LSST', 'HUBBLE'], normed=True)

  for i, band in enumerate(filters.keys()):
    pl.plot(filters[band]['ls'], filters[band]['Ts'], label=band)

  pl.legend(ncol=4)
  pl.show()

  print("\n\nDone.\n\n")
