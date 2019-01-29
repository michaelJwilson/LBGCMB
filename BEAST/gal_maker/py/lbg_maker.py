import  os
import  numpy                  as       np
import  pylab                  as       pl
import  astropy.constants      as       const
import  astropy.units          as       u 

from    astropy.io             import   fits  
from    utils                  import   latexify
from    astropy                import   units               as  u
from    redshift_spectra       import   redshift_spectra
from    prep_filters           import   prep_filters
from    scipy.interpolate      import   interp1d
from    app_mags               import   get_appmags
from    astropy.table          import   Table
from    desispec.interpolation import   resample_flux
from    astropy.cosmology      import   FlatLambdaCDM
from    collections            import   OrderedDict
from    uv_filter              import   get_detband


cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

def lbg_maker(ngal=None, restframe=False, printit=False, test=False, seed=314, redshifts=None, magnitudes=None, hiwave=1.4e4):
   np.random.seed(seed)

   if type(redshifts)  == type(None):
     redshifts   =   1.5 + np.linspace(0.,   4.0,    20)

   if type(magnitudes) == type(None):
     magnitudes  =  20.0 +   np.arange(0.,   7.0,   0.1)

   if restframe:
     ##  EW magnitudes for quantiles of Shapley spectra;  Left col. of pg. 8 (https://arxiv.org/pdf/astro-ph/0301230.pdf).  Includes absorption edge. 
     eqv_width     = [       -14.92,         -1.10,         11.00,         14.30,         52.63]
     spectra       = ['Q0_199.fits', 'Q2_198.fits', 'Q3_199.fits', 'QX_811.fits', 'Q4_198.fits']
   
   else:
     ##  No need to sample both Q3 and QX in LBG generation. 
     eqv_width     = [       -14.92,         -1.10,         11.00,         52.63]
     spectra       = ['Q0_199.fits', 'Q2_198.fits', 'Q3_199.fits', 'Q4_198.fits']

   ##  Order by decreasing EW.
   eqv_width     =  eqv_width[::-1]
   spectra       =    spectra[::-1]

   eqv_width     =  OrderedDict(zip(spectra, eqv_width))

   ##  Cut on EW considered.
   ##  spectra   =   spectra[:2]
   
   if printit:
     print('\n\nSolving for redshifts:\n'  + ''.join('%.3lf, ' % x for x in redshifts))
     print('\n\nSolving for magnitudes:\n' + ''.join('%.3lf, ' % x for x in magnitudes))
   
   root        =  os.environ['BEAST'] + '/gal_maker/dat/composites/'
   
   filters     =  prep_filters()
   
   ##  Used to set detection band. 
   lsstfilters =  prep_filters(['LSST'])
   
   wave        =  np.arange(3.6e3, hiwave, 0.2)  ##  Matched to gal. maker. 
   count       =  0

   flux        =  []
   
   ##  Meta.
   objtypes    =  [] 
   subtypes    =  []
   meta_zs     =  []
   meta_eqws   =  []
   meta_mags   =  []
   meta_dbands =  []

   nmax        =  len(redshifts) * len(spectra) * len(magnitudes)

   if (ngal is None):
     ngal      =  nmax

   elif(ngal > nmax):
     ngal      =  nmax

     raise ValueWarning('Requested LBG number is larger than that available.')

   else:
     pass

   runs = [[z, spec, mag] for mag in magnitudes for z in redshifts for spec in spectra]
   
   if test:
     ##  EW magnitudes for quantiles of Shapley spectra;  Left col. of pg. 8 (https://arxiv.org/pdf/astro-ph/0301230.pdf).  Includes absorption edge.                                                                                                                                      
     eqv_width     =  [       -14.92,         -1.10,         11.00,         14.30,         52.63]
     spectra       =  ['Q0_199.fits', 'Q2_198.fits', 'Q3_199.fits', 'QX_811.fits', 'Q4_198.fits'] 
     eqv_width     =  OrderedDict(zip(spectra, eqv_width))

     runs          =  [[5.5, spectra[0], 23.5], [4.5, spectra[2], 23.5], [3.5, spectra[3], 23.5], [2.5, spectra[4], 23.5]]
     printit       =  True
     ngal          =  4

   if restframe:
     ##  One for each Shapley composite.                                                                                                                 
     ngal       =  5
     runs       =  [[0.0, spec, 23.5] for spec in spectra]
     printit    =  True

     print(runs)

   ##  Cut on number of galaxies requested.
   runs = np.array(runs)
   runs = runs[np.random.choice(runs.shape[0], ngal, replace=False)]
   
   ## flux_factors = 10.0 ** (0.4 * (magnitudes[0] - magnitudes))
   for run in runs:
     redshift  = np.float(run[0])
     specname  =          run[1]
     magnitude = np.float(run[2])

     if printit:
       print('Solving for redshift: %.3lf, mag: %.3lf, spec: %s' % (redshift, magnitude, specname)) 
     
     fname = root + specname
     
     if fname.split('.')[-1] == 'fits': 
       dat     =   fits.open(fname)

       ##  WAT1_001 = 'wtype=linear label=Wavelength units=Angstroms'.  See http://www-wfau.roe.ac.uk/6dFGS/fits.html
       waves   =   u.AA * (dat[0].header['CRVAL1'] - (dat[0].header['CRPIX1'] - np.arange(dat[0].header['NAXIS1'])) * dat[0].header['CDELT1'])

       vs      =   waves.to(u.Hz, equivalencies = u.spectral())         ##  AA to Hz.
       Fv      =   dat[0].data * u.erg / u.s / u.cm / u.cm / u.Hz       ##  [erg / s / sq. cm / Hz].

       ##  Flux enforced to be positive definite. 
       Fv[Fv < 0.0] = 0.0 * u.erg / u.s / u.cm / u.cm / u.Hz

     else:
       dat     =   np.loadtxt(fname)

       waves   =   u.AA * dat[0]
       vs      =   waves.to(u.Hz, equivalencies = u.spectral())         ##  AA to Hz.

       Fl      =   dat[1] * 1.e-17 * u.erg / u.s / u.cm / u.cm / u.AA
       Fv      =   Fl.to(u.erg / u.s / u.cm / u.cm / u.Hz, equivalencies = u.spectral_density(waves))

     ##  Extend break and flat Fv.
     redend  =  np.mean(Fv.value[(1.6e3 < waves.value) & (waves.value < 1.8e3)])

     dwave   =  waves[1].value - waves[0].value

     lowext  =  np.arange(1.e2, waves.value.min() + dwave, dwave)
     hiext   =  np.arange(waves.value.max() + dwave, 2.e4, dwave)

     waves   =  np.concatenate([lowext, waves.value, hiext]) * u.AA
     Fv      =  np.concatenate([np.zeros_like(lowext), Fv.value, redend * np.ones_like(hiext)]) * u.erg / u.s / u.cm / u.cm / u.Hz
     vs      =  waves.to(u.Hz, equivalencies = u.spectral())

     if not restframe:
       ## Redshift spectra for quickspectra input.  'Else', leave unredshifted for redrock template.
       vs, Fv  =   redshift_spectra(vs.value, Fv.value, redshift)

       vs      =   vs * u.Hz
       Fv      =   Fv * u.erg / u.s / u.cm / u.cm / u.Hz

       waves   =   vs.to(u.AA, equivalencies = u.spectral())            ##  Redshifted wavelength.
 
     print(vs.value)
       
     exit(0)

     mags  =  get_appmags(vs.value, Fv.value, filters, printit = False)
     dband =  get_detband(redshift, lsstfilters, printit=False, restframe=restframe)
     
     Fv   *=  10. ** (0.4 * (mags[dband] - magnitude))                  ##  Normalise to magnitudes[0], e.g. 20.0 in g.

     ##  And check.
     ##  mags  = get_appmags(vs, Fv, filters, printit = True)
       
     ##  F_lambda from F_nu.
     Fl   = Fv.to(u.erg / u.s / u.cm / u.cm / u.AA, equivalencies = u.spectral_density(waves))       ##  [       ergs / s / sq. cm / A]
       
     if restframe:
       wave  = waves.value
       Fl   /= 1.e-17  

     else:
       Fl   /= 1.e-17                                                                                ##  [1.e-17 ergs / s / sq. cm / A]        
       Fl    = resample_flux(wave, waves.value, Fl.value) * u.erg / u.s / u.cm / u.cm / u.AA

     ##  Magnitude is arbitrary for redrock template. 
     objtypes.append('LBG')
     subtypes.append(specname.split('_')[0])
     meta_zs.append(redshift)
     meta_eqws.append(eqv_width[specname])
     meta_mags.append(magnitude)
     meta_dbands.append(dband)
        
     flux.append(Fl.value)

     count += 1

     if count >= ngal:
       flux =  tuple(flux)
       flux =  np.column_stack(flux)
       flux =  flux.T

       meta =  Table([objtypes, subtypes, meta_zs, meta_eqws, meta_mags, meta_dbands], names=('OBJTYPE', 'SUBTYPE', 'REDSHIFT', 'Lya-EW', 'AB mag', 'dbands'),\
                      meta={'name': 'LBG MAKER META'})

       if flux.min() < 0.0:
         raise ValueError('ERROR:  z: %.3lf, spectra: %s' % (redshift, specname))

       return  flux, wave, meta


if __name__ == '__main__':
  print('\n\nWelcome to (Shapley) LBG maker.')

  ngal              =      15
  test              =   False
  target_type       =    'lbg'
  restframe         =    True
  linecalc          =   False
  save              =   False

  meta_name         = os.environ['BEAST']   + '/gal_maker/dat/Tables/galmaker-%s-meta-nir.txt' % target_type
  infile            = os.environ['SCRATCH'] + '/desi/simspec/%s-input-spectra-nir.fits' % target_type

  redshifts         =   1.5 + np.linspace(0.,   4.0,    20)
  magnitudes        =  20.0 +   np.arange(0.,   7.0,   0.1)
  
  flux, wave, meta  =   lbg_maker(ngal=ngal, restframe=restframe, printit=False, test=test, redshifts=redshifts, magnitudes=magnitudes)
  
  if test:
    for i, x in enumerate(flux):
      pl.semilogy(wave, x, label=str(meta['r'][i]) + ', %.1lf' % meta['REDSHIFT'][i])

    pl.ylim(0.1, 10.1)
    pl.legend(ncol=2, loc=1)
    pl.savefig('../plots/lbg_spectra.pdf')
  
  print(wave.shape)
  print(flux.shape)
  print(meta)

  dwave = wave[1] - wave[0]

  if linecalc:
    ## ---  Ly-alpha line calc.  --- ##
    for i, xx in enumerate(flux[(meta['SUBTYPE'] == 'Q4') & (meta['REDSHIFT'] > 2.1) & (meta['REDSHIFT'] < 2.2)]):
       pl.clf()

       redshift      = meta['REDSHIFT'][i]
       lumdist       = cosmo.luminosity_distance(z=meta['REDSHIFT'][i]).to(u.cm).value

       Lya           = 1215.24 * (1. + redshift)
       dlambda       =     5.8 * (1. + redshift)
    
       linecut       = ((Lya - dlambda) <= wave) & (wave <= (Lya + dlambda))
       ## continuumcut = ((Lya + dlambda) <= wave) & (wave <= (Lya + dlambda + 1.))

       wavecut       = wave[linecut]
       fluxcut       =   xx[linecut]
       fluxcut      *= 1.e-17
       
       lineflux      = np.sum(fluxcut) * dwave
       contflux      = np.sum(np.ones_like(fluxcut)) * np.mean([fluxcut[0], fluxcut[-1]]) * dwave

       contcorrected = lineflux - contflux

       print('%.6le' % (np.sum(fluxcut) * dwave))

       pl.axvline(x=Lya,     ymin=0., ymax=1., c='k')
       pl.axhline(y=np.mean([fluxcut[0], fluxcut[-1]]), xmin=0., xmax=1., c='k')
       
       pl.plot(wave,  xx * 1.e-17, 'k-', alpha=0.4)
       pl.plot(wavecut,   fluxcut,  'o', markersize=4, label=str(i))
       pl.title('%.3le, %.3le [ergs/s/cm/cm], %.3lf [ergs/s]' % (lineflux, contcorrected, np.log10(contcorrected * lumdist * lumdist)))
       
       pl.xlim([Lya - 25., Lya + 25.])
       pl.legend()
       pl.show()

       break

  if save:
    meta.write(os.environ['BEAST'] + '/gal_maker/dat/Tables/galmaker-%s-meta.txt' % target_type, format='ascii', overwrite=True)

    hdr    = fits.Header()

    hdr['EXTNAME'] = 'WAVELENGTH'
    hdr['BUNIT']   = 'Angstrom'

    fits.writeto(infile, wave, header=hdr, overwrite=True)

    hdr['EXTNAME'] = 'FLUX'
    hdr['BUNIT']   = '10^-17 erg/(s*cm^2*Angstrom)'  # Satisifes FITS standard AND Astropy-compatible.

    fits.append(infile, flux, header=hdr)
  
    hdr['EXTNAME'] = 'REDSHIFTS'
    hdr['BUNIT']   = 'DIMENSIONLESS'                 # Satisifes FITS standard AND Astropy-compatible.                                                                                                                             
    fits.append(infile, redshifts, header=hdr)

    hdr['EXTNAME'] = 'MAGNITUDES'
    hdr['BUNIT']   = 'DIMENSIONLESS'                 # Satisifes FITS standard AND Astropy-compatible.                                                                                                                           

    fits.append(infile, magnitudes, header=hdr)
    

  print('\n\nDone.\n\n')
