import  collections

import  os
import  numpy             as np
import  astropy.constants as const


def prep_filters(names=['LSST', 'VIDEO'], normed=False):
  root     = os.environ['LBGCMB'] + '/dropouts/filters/'
  names    = [x.upper() for x in names]
  filters  = collections.OrderedDict()                                              ## Generate the required filters; Euclid: Y, J and H.                    

  if 'LSST' in names:
    for band in ['u', 'g', 'r', 'i', 'z', 'y']:
      ppkey  = r"LSST-$" + "%s" % band + r"$"                                       ## Pretty print version of key.  

      fname  = root + 'lsst/total_%s.dat' % band

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
      for band in ['Y', 'J', 'H', 'K']:
        ppkey  = r"VIDEO-$" + "%s" % band + r"$"                                       ## Pretty print version of key.

        fname  = root + "video/%s.txt" % band
        data   = np.loadtxt(fname)

        ls     = data[:,0]
        vs     = (1e10/ls) * const.c.to('m/s').value

        if normed:
          filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1] / data[:,1].max()}

        else:
          filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1]}

      names.remove('VIDEO')

  if 'SUBARU' in names:
    print('Adding Subaru filters (B and V).')

    for band in ['B', 'V']:
      ppkey  = r"SUBARU-$" + "%s" % band + r"$"
      fname  = root + "subaru/%s.pb" % band
      data   = np.loadtxt(fname)

      ls     = data[:,0]
      vs     = (1.e10 / ls) * const.c.to('m/s').value

      if normed:
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1] / data[:,1].max()}

      else:
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1]}

    names.remove('SUBARU')

  if 'EUCLID' in names:
    ##  https://arxiv.org/pdf/1710.09585.pdf
    print('Adding Euclid filters (VIS (riz), Y, J, H).')

    for band in ['VIS', 'Y', 'J', 'H']:
      ppkey  = r"EUCLID-$" + "%s" % band + r"$"
      fname  = root + "euclid/%s.pb" % band
      data   = np.loadtxt(fname)

      ls     = data[:,0] ## Angstroms. 
      vs     = (1.e10 / ls) * const.c.to('m/s').value

      if normed:
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1] / data[:,1].max()}

      else:
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1]}

    names.remove('EUCLID')

  if 'JKC' in names:
    print('Adding JKC filters (I).')

    for band in ['I']:
      ppkey  = r"JKC-$" + "%s" % band + r"$"
      fname  = root + "jkc/%s.pb" % band
      data   = np.loadtxt(fname)

      ls     = data[:,0]
      vs     = (1.e10 / ls) * const.c.to('m/s').value

      if normed:
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1] / data[:,1].max()}

      else:
        filters[band] = {'ppkey': ppkey, 'fname': fname, 'ls': ls, 'vs': vs, 'Ts': data[:,1]}

    names.remove('JKC')

  if len(names) != 0:
    ##  Catch any requested filters that are not either LSST or VIDEO.                                                                                                                                                                     
    raise ValueWarning('Requested filters are not available:' + '  '.join(x for x in name))
  
  return  filters
  

if __name__ == "__main__":
  import pylab as pl


  print('\n\nWelcome to prep-filters.\n\n')

  filters = prep_filters(['Euclid'])

  for i, band in enumerate(filters.keys()):
    print(band)

    pl.plot(filters[band]['ls'], filters[band]['Ts'], label = filters[band]['ppkey'])

  ##  for filter in filters.values():
  ##    print filter['ppkey']

  

  pl.legend(loc=2, ncol=4, frameon=False)

  pl.ylim(0.0, 1.2)
  pl.savefig('plots/filters.pdf')

  print("\n\nDone.\n\n")
