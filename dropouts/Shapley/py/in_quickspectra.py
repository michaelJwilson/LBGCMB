import  os
import  numpy                as       np
import  pylab                as       pl
import  astropy.constants    as       const

from    astropy.io           import   fits  
from    utils                import   latexify
from    astropy              import   units               as  u
from    redshift_spectra     import   redshift_spectra
from    prep_filters         import   prep_filters
from    scipy.interpolate    import   interp1d
from    app_mags             import   get_appmags
from    rrtemplate_io        import   create_template


def get_filters(plotit=True):
  filters = prep_filters()

  if plotit:
    for i, band in enumerate(['r']):
      pl.fill_between(filters[band]['ls'], 5. * filters[band]['Ts'], alpha=0.3, label=filters[band]['ppkey'], color='g')

  return  filters

def get_noisecurve(exp, plotit = True):
  fname  =  os.environ['LBGCMB'] + '/quickspectra/dat/quickspectra/nosource/real0/nosource_exp' + str(exp) + '.fits'
  noise  =  fits.open(fname)

  pl.plot(noise["B_WAVELENGTH"].data, 1. / noise["B_IVAR"].data[0], label='', c='k', alpha=0.9)
  pl.plot(noise["R_WAVELENGTH"].data, 1. / noise["R_IVAR"].data[0], label='', c='k', alpha=0.9)
  pl.plot(noise["Z_WAVELENGTH"].data, 1. / noise["Z_IVAR"].data[0], label='', c='k', alpha=0.9)

  pl.title('%.0lf mins.' % (exp / 60.))

  return  noise


if __name__ == '__main__':
   print('\n\nWelcome to a constructor for a redrock template and\nquickspectra input from Shapley composite.\n\n')

   ##  EW magnitudes for quantiles of Shapley spectra;  Left col. of pg. 8 (https://arxiv.org/pdf/astro-ph/0301230.pdf).  Includes absorption edge. 
   eqv_width   = [       -14.92,         -1.10,         11.00,         14.30,         52.63]
   spectra     = ['Q0_199.fits', 'Q2_198.fits', 'Q3_199.fits', 'QX_811.fits', 'Q4_198.fits']

   ##  Order by decreasing EW.
   eqv_width   =  eqv_width[::-1]
   spectra     =    spectra[::-1]

   sample_rate =  1

   redshifts   =   1.5 + np.linspace(0., 1.50, 15)[::sample_rate]
   magnitudes  =  20.0 +   np.arange(0., 6.00,  1)

   print('\n\nSolving for redshifts:\n'  + ''.join('%.3lf, ' % x for x in redshifts))
   print('\n\nSolving for magnitudes:\n' + ''.join('%.3lf, ' % x for x in magnitudes))

   ##  For redrock template creation, redshifted should be False.
   redshifted  =    True
   plotit      =    False

   fnames      =  [os.environ['LBGCMB'] + '/dropouts/Shapley/dat/composites/' + file for file in spectra]
   fnames     +=  [os.environ['LBGCMB'] + '/quickspectra/spectra/spec-qso-z1.5-rmag22.24.dat']
   fnames     +=  [os.environ['LBGCMB'] + '/quickspectra/spectra/spec-elg-o2flux-8e-17.dat']

   latexify(columns=2, fontsize=10, ggplot=True)

   ##  Get noise curves and filters. 
   ##  exp     =  60 * 60                             ##  [Seconds].
   ##  noise   =  get_noisecurve(exp, plotit=plotit)

   filters     =  get_filters(plotit=plotit)

   for redshift in redshifts:
    print('\n\nSolving for redshift: %.3lf.' % redshift) 

    owaves     =   np.arange(1., 1.e4, 0.1)
    output     =  [owaves]

    for kk, fname in enumerate(fnames):
     print('\nLoading:  %s' % fname)

     if fname.split('.')[-1] == 'fits':
       dat     =   fits.open(fname)

       ##  WAT1_001 = 'wtype=linear label=Wavelength units=Angstroms'.  See http://www-wfau.roe.ac.uk/6dFGS/fits.html
       waves   =   u.AA * (dat[0].header['CRVAL1'] - (dat[0].header['CRPIX1'] - np.arange(dat[0].header['NAXIS1'])) * dat[0].header['CDELT1'])

       vs      =   waves.to(u.Hz, equivalencies = u.spectral())         ##  AA to Hz.
       Fv      =   dat[0].data * u.erg / u.s / u.cm / u.cm / u.Hz       ##  [erg / s / sq. cm / Hz].

     else:
       dat     =   np.loadtxt(fname)

       waves   =   u.AA * dat[0]
       vs      =   waves.to(u.Hz, equivalencies = u.spectral())         ##  AA to Hz.

       Fl      =   dat[1] * 1.e-17 * u.erg / u.s / u.cm / u.cm / u.AA
       Fv      =   Fl.to(u.erg / u.s / u.cm / u.cm / u.Hz, equivalencies = u.spectral_density(waves))

     ##  Create interpolator and extrapolate flat Fv to 1 um with fill_value.                                                                         
     Fv        =   interp1d(vs.value, Fv.value, kind='linear', bounds_error=False, fill_value=Fv[-1], assume_sorted=False)

     waves     =   np.arange(np.min(waves.value), 1.5e4, 0.1) * u.AA    ##  [Angstroms].                                                       
     vs        =   waves.to(u.Hz, equivalencies = u.spectral())

     Fv        =   Fv(vs.value) * u.erg / u.s / u.cm / u.cm / u.Hz      ##  [erg / s / sq. cm / Hz].

     if redshifted:
       ## Redshift spectra for quickspectra input.  'Else', leave unredshifted for redrock template.
       vs, Fv  =   redshift_spectra(vs, Fv, redshift)
       waves   =   vs.to(u.AA, equivalencies = u.spectral())            ##  Redshifted wavelength.

       mags    =   get_appmags(vs, Fv, filters, printit = False)
       Fv     *=   10. ** (0.4 * (mags['r'] - magnitudes[0]))           ##  Normalise to 21.5 mag. in g.

     ##  And check.
     ##  mags  =   get_appmags(vs, Fv, filters, printit = True)
    
     ##  F_lambda from F_nu.
     Fl     = Fv.to(u.erg / u.s / u.cm / u.cm / u.AA, equivalencies = u.spectral_density(waves))    ##  [       ergs / s / sq. cm / A]
     Fl    /= 1.e-17                                                                                ##  [1.e-17 ergs / s / sq. cm / A] 
        
     Fl     = interp1d(waves.value, Fl.value, kind='linear', bounds_error=False, fill_value=0.0, assume_sorted=False)
     Fl     = Fl(owaves)

     if redshifted:
       flux_factors  =   10. ** (0.4 * (magnitudes[0] - magnitudes))
       output       +=  [Fl * fluxfactor for fluxfactor in flux_factors]
     
     else:
       ##  Twinned with unredshifted spectra for redrock template; .fits output.                                                                     
       create_template(owaves, Fl, file.split('.')[0], printit = False)

       ##  Only create one (rest-frame) template.
       exit(0)

     if plotit:
         pl.plot(owaves + 50. * kk, Fl, label='% .1lf' % eqv_width[kk] + r'$\AA$')

         pl.xlabel(r'Wavelength [$\AA$]')
         pl.ylabel(r'F$_\lambda$ [$10^{-17}$ ergs/$s$/cm$^2$/$\AA$]')

         pl.xlim(0.e0, 1.e4)
         pl.legend(ncol=1)

         pl.savefig(os.environ['LBGCMB'] + '/dropouts/Shapley/plots/composites.pdf', bbox_inches='tight')

    output =  tuple(output)
    output =  np.column_stack(output)

    ##  Twinned with redshifted spectra for quickspectra input; .txt output.                                                                                
    np.savetxt(os.environ['LBGCMB'] + '/quickspectra/dat/in_quickspectra/baseline/spec-all-z%.3lf.dat' % (redshift), output, fmt='%.6le')

    ##  np.savetxt(os.environ['LBGCMB'] + '/quickspectra/dat/in_quickspectra/Shapley/spec-Shapley-%s-z%.3lf.dat' % (file.split('.')[0], redshift),\
    ##                                                                                                              output, fmt='%.6le')
    
print('\n\nDone.\n\n')
