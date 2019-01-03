
import  re
import  os 
import  numpy            as      np
import  astropy.io.fits  as      fits

from    glob             import  glob
from    utils            import  latexify


def native_endian(data):
    '''
    Convert numpy array data to native endianness if needed.                                                                                           
    Returns new array if endianness is swapped, otherwise returns input data.
    
    Context:
    By default, FITS data from astropy.io.fits.getdata() are not Intel
    native endianness and scipy 0.14 sparse matrices have a bug with
    non-native endian data.

    ##  https://github.com/desihub/desispec io.util.native_endian
    '''

    if data.dtype.isnative:
      return data

    else:
      return data.byteswap().newbyteorder()

def create_template(rwavelength, flux, type, printit = False, root=None):
    '''                                                                                                                                                     
    Given type, linearly sampled restframe wavelength [Angstroms] and 
    corresponding restframe flux [1e-17 ergs/s/cm^2/A], write a rrtemplate-type.fits 
    file.                                                                                          
    '''
    if not np.all(np.diff(rwavelength) > 0):
      raise ValueError('\n\nProvided wavelength is not monotonically increasing.') 

    hdr             = fits.Header()
    
    hdr['LOGLAM']   = 0                                ## Linear spacing. 
    hdr['CRVAL1']   = rwavelength[0]                   ## Starting wavelength [angstroms]
    hdr['CDELT1']   = rwavelength[1] - rwavelength[0]  ## dlambda; assumes linear.
    hdr['RRTYPE']   = type.upper()

    hdr['EXTNAME']  = 'BASIS_VECTORS'

    if type == 'LBG':
        hdr['RRSUBTYP'] = '' ## 'NoLines'

    else:
        hdr['RRSUBTYP'] = ''

    if root is None:
        root = os.environ['LBGCMB'] + '/redrock/templates/ext_wave/'

    flux            = flux.reshape(1, len(flux))
    flux            = native_endian(flux)

    hdu             = fits.PrimaryHDU(flux, header=hdr)

    hdu.writeto(root + 'rrtemplate-%s.fits' % type, overwrite = True)

    print("Template written.")

    if printit is True:
      print('LOGLAM:  %s' % str(hdr['LOGLAM']))  
      print('CRVAL1:  %s' % str(hdr['CRVAL1']))
      print('CDELT1:  %s' % str(hdr['CDELT1']))

      for i, wave in enumerate(rwavelength[::10]):
        print "%.6lf \t %.6le" % (wave, flux[0][10 * i])

def read_template(type='qso', printit=False, root=None):
    '''
    Read a redrock template and return the restframe wavelength and 
    flux [1e-17 ergs/s/A/cm^2]    
    '''

    import os 

    if root is None:
        root = os.environ['LBGCMB'] + '/redrock/templates/ext_wave/'

    fname    = root + "rrtemplate-%s.fits" % type

    print("\n\nReading redrock template:  %s\n\n" % fname)

    template = fits.open(fname)

    if printit:
        print  template.info()        
        print  template[0].header  
    
    ntemp    = template[0].data.shape[0]
    nwave    = template[0].data.shape[1]

    hdr      = template[0].header

    dlambda  = hdr['CDELT1']
    length   = hdr['NAXIS1']
    lmin     = hdr['CRVAL1']

    wave     = lmin + dlambda * np.arange(length)

    if 'LOGLAM' in hdr and hdr['LOGLAM'] != 0:
      wave   = 10.**wave

    result   = {}

    result['ntemp'] = ntemp
    result['wave']  = wave
    
    for i in np.arange(ntemp):
        result['temp_%d' % i] = template[0].data[i, :]
        
    return  result


if __name__ == '__main__':
    import  pylab              as      pl
    import  matplotlib.pyplot  as      plt

    from    pylab              import  rcParams


    print("\n\nWelcome to rrtemplate_io.")

    ## data     = np.loadtxt(os.environ['LBGCMB'] + '/quickspectra/spectra/BC03/restframe/spec-BC03-z0.0.dat')
    ## create_template(data[:,0], data[:,1], type = 'LBG', printit = False)

    root     = None
    root     = os.environ['LBGCMB'] + '/redrock/templates/ext_wave/'
    
    types    = {'Q2_198': None, 'LBG': None, 'ELG': None}
    ## types = {'galaxy': None,  'qso':    None, 'LBG': None, 'star-B': None}
    ## types = {'star-M': None,  'star-K': None, 'star-G': None, 'star-F': None,  'star-A': None, 'star-B': None}

    latexify(fig_width=12., fig_height=6., equal=False, fontsize=10, ggplot=True, usetex=True)
    
    pl.axvline(x=1216., ymin=0., ymax=1., c='k')

    for type in types:
        types[type] = read_template(type=type, printit=True, root=root)
        
        for i in np.arange(types[type]['ntemp']):
          index     = np.where(np.abs(types[type]['wave'] - 6.e3) == np.abs(types[type]['wave'] - 6.e3).min())     
          norm      = types[type]['temp_%d' % i][index]

          pl.plot(types[type]['wave'], types[type]['temp_%d' % i] / norm, lw=.5, label='Template %s:  %d' % (re.escape(type), i))

        pl.xlim(types[type]['wave'].min(), types[type]['wave'].max())
        ## pl.ylim(-50., 100.)

        pl.axhline(y=0., xmin=0., xmax=1., c='y', lw=0.2)

        pl.xlabel(r'Wavelength $[\AA]$')
        pl.ylabel(r'$F_{\lambda} \ [10^{-17} \ \rm{ergs}/s/\AA/\rm{cm}^2$]')

        pl.xlim(1.e0, 1.5e4)

        pl.title(type.upper())

        pl.legend(loc=1, ncol=4)

    pl.savefig("plots/rrtemplate.pdf", bbox_inches = "tight")
         
    print("\n\nDone.\n\n")
