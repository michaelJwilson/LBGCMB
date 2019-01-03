import  os
import  numpy              as      np 
import  pylab              as      pl
import  astropy.io.fits    as      fits

from    rrtemplate_io      import  create_template
from    scipy.interpolate  import  interp1d 


## Stellar type grouping by temperature; templates do not cover o-type with Teff > 30k.                                                                     
stellar_classes  = {'B': [1.0e4, 3.0e4], 'A': [7.5e3, 1.0e4], 'F': [6.0e3,  7.5e3],
                    'G': [5.2e3, 6.0e3], 'K': [3.7e3, 5.2e3], 'M': [2.4e3,  3.7e3]}

if __name__ == "__main__":
    print("\n\nWelcome to a conversion of basis_templates to redrock templates.\n\n")

    root      = os.environ['DESI_BASIS_TEMPLATES']

    '''
    type, ver =  'elg', 2.0
    type, ver =  'lrg', 2.0
    type, ver =  'qso', 2.0
    '''
    type, ver = 'star', 2.2

    fpath     = root + type + '_templates_v%.1lf.fits' % ver
    data      = fits.open(fpath) 

    print  data.info()

    for subtype in stellar_classes:
      lo_temp   = stellar_classes[subtype][0]
      hi_temp   = stellar_classes[subtype][1]

      temps     = data[1].data['TEFF']
      cut       = (lo_temp < temps) & (temps < hi_temp)

      wave      = data[2].data

      flux      = data[0].data[cut, :] 
      flux      = flux[0, :]

      flux_int  = interp1d(wave, flux, kind='linear', fill_value=0.0, assume_sorted=False, bounds_error=None)

      wave      = np.arange(wave.min(), 6.e4, 1.)
      flux      = flux_int(wave)

      ## Create redrock template extended in wavelength. 
      create_template(wave, flux, type = type.upper() + '-' + subtype, printit = False, root='./redrock/')

    print("\n\nDone.\n\n")
