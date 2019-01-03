import  os
import  numpy                  as       np
import  pylab                  as       pl
import  astropy.constants      as       const

from    astropy.io             import   fits  
from    utils                  import   latexify
from    astropy                import   units             as  u
from    redshift_spectra       import   redshift_spectra
from    scipy.interpolate      import   interp1d
from    scipy.ndimage.filters  import   gaussian_filter


latexify(columns=2, fontsize=10, ggplot=True)

##  Plot noise curves for given exposure lengths. 
for ii, exp in enumerate(60 * 60 * np.arange(1, 6, 1)):
  for channel, color in zip(['B', 'R', 'Z'], ['b', 'g', 'r']):
    ##  Plot (pure) noise curves.
    fname       =   os.environ['LBGCMB'] + '/quickspectra/dat/nosource/real0/nosource_exp' + str(exp) + '.fits'
    noise       =   fits.open(fname)

    wave        =   noise[channel + '_WAVELENGTH'].data
    smooth      =   gaussian_filter(1. / noise[channel + '_IVAR'].data[0], sigma=3)

    pl.plot(wave[::5], smooth[::5], label='', c=color, alpha=1.0 - 0.1 * ii)

    print exp / 60. / 60.

pl.xlabel(r'Wavelength [$\AA$]')
pl.ylabel(r'N$_\lambda$ [$10^{-17}$ ergs/$s$/cm$^2$/$\AA$]') 

pl.yscale('log')

pl.xlim(3.6e3, 1.e4)
pl.legend(ncol=1)

pl.savefig('plots/exposures.pdf', bbox_inches='tight')
