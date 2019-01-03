from    __future__            import  division

import  astropy.io.fits       as      pyfits
import  pylab                 as      pl
import  numpy                 as      np
import  matplotlib.pyplot     as      plt
import  os

from    matplotlib.offsetbox  import  TextArea, AnnotationBbox
from    utils                 import  latexify


def join_arms(exp=None, plotit=False):
    exp                  = '%.0lf' % exp
    
    fname                = 'dat/nosource/real2/nosource_exp' + exp + '.fits'

    print('\n\nLoading %s.' % fname)

    exposure             = np.float(exp)
    exposure            /= 60.  ## Minutes                                                                                                               

    title                = fname.split('/')[-1]

    ## [type, exp]       = title.split('-')

    ## [spec, type, z, exp] = title.split('-')                                                                                                         
    ## [spec, type, z, mag] = title.split('-')                                                                                                           
    h                    = pyfits.open(fname)

    h.info()
    
    brjoin               = 0.5 * (h["B_WAVELENGTH"].data.max() + h["R_WAVELENGTH"].data.min())
    rzjoin               = 0.5 * (h["R_WAVELENGTH"].data.max() + h["Z_WAVELENGTH"].data.min())

    bflux                = h["B_FLUX"].data[0][h["B_WAVELENGTH"].data < brjoin]
    rflux                = h["R_FLUX"].data[0][h["R_WAVELENGTH"].data > brjoin]
    
    bwave                = h["B_WAVELENGTH"].data[h["B_WAVELENGTH"].data < brjoin]
    rwave                = h["R_WAVELENGTH"].data[h["R_WAVELENGTH"].data > brjoin]
    
    rflux                = rflux[rwave < rzjoin]
    zflux                = h["Z_FLUX"].data[0][h["Z_WAVELENGTH"].data > rzjoin]
    
    rwave                = rwave[rwave < rzjoin]
    zwave                = h["Z_WAVELENGTH"].data[h["Z_WAVELENGTH"].data > rzjoin]

    waves                = np.concatenate([bwave, rwave, zwave])
    flux                 = np.concatenate([bflux, rflux, zflux])

    ## pl.savetxt(os.environ['FILE'] + '-exp' + os.environ['EXPT'] + '-singlearm' + '.txt', np.column_stack((waves, flux)), fmt=%.6le)

    if plotit:
      latexify(columns=2, equal=False, fontsize=10, ggplot=True)
        
      ## pl.plot(waves, flux, label='%.1lf mins.' % exposure)

      pl.xlim(3.5e3, 1.e4)
      pl.ylim(0.,      6.)

      pl.xlabel(r"$\lambda$ [$\AA$]")
      pl.ylabel(r"$F_\lambda$ / $10^{-17}$ / ergs / $s$ / cm$^2$ / $\AA$")

      pl.legend()

      pl.plot(h["B_WAVELENGTH"].data, 1. / h["B_IVAR"].data[0], label='%.1lf mins.' % exposure)
      pl.plot(h["R_WAVELENGTH"].data, 1. / h["R_IVAR"].data[0], label='')
      pl.plot(h["Z_WAVELENGTH"].data, 1. / h["Z_IVAR"].data[0], label='')

      pl.savefig('plots/nosource/ivar_exp_%.0lf.pdf' % exposure)
    

if __name__ == "__main__":
    print('\n\nWelcome to join arms.\n\n')
    
    exposures  = 60. * 15. * np.arange(1, 22, 1)     ## Seconds.  

    for exp in exposures[::5]:
      join_arms(exp=exp, plotit=True)
    
    print('\n\nDone.\n\n')
