import  os
import  glob
import  numpy                  as       np
import  pylab                  as       pl
import  astropy.constants      as       const
import  astropy.units          as       u 
import  matplotlib.pyplot      as       plt

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
from    read_ised              import   read_ised
from    numpy.random           import   choice
from    extinction             import   calzetti00
from    extinction             import   apply               as ext_apply
from    madau                  import   lephare_madau
from    uv_filter              import   get_detband


root  = os.environ['BEAST']
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

def BC03_maker(ngal=None, restframe=False, printit=False, test=False, seed=314, redshifts=None, magnitudes=None, alliseds=False, calzetti=False, madau=False):
   np.random.seed(seed)
   
   if type(redshifts)  == type(None):
     redshifts   =   1.5 + np.linspace(0.,   4.0,    20)

   if type(magnitudes) == type(None):
     magnitudes  =  20.0 +   np.arange(0.,   7.0,   0.1)
   
   if printit:
     print('\n\nSolving for redshifts:\n'  + ''.join('%.3lf, ' % x for x in redshifts))
     print('\n\nSolving for magnitudes:\n' + ''.join('%.3lf, ' % x for x in magnitudes))

     
   filters     =  prep_filters()
   lsstfilters =  prep_filters(['LSST'])

   ##  wave    =  np.arange(1.0e0, 1.e4, 0.1) 
   owave       =  np.arange(3.6e3, 1.e4, 0.1)  ##  Matched to gal. maker. 
   count       =  0

   flux        =  []
   
   ##  Meta.
   objtypes    =  [] 
   subtypes    =  []
   meta_zs     =  []
   meta_dbands =  []
   meta_ages   =  []
   meta_imfs   =  []
   meta_taus   =  []
   meta_metals =  []
   meta_mags   =  []
   meta_calz   =  []

   root        =  os.environ['BEAST'] + '/gal_maker/dat/BC03/salpeter/'
   spectra     =  []

   ##  Save metallicity, tau and ages for every retained template. 
   Ms          =  []
   Is          =  []
   Ts          =  []
   Ages        =  []


   files       =  []

   if alliseds:
      files    =  glob.glob(os.environ['BEAST'] + '/gal_maker/dat/BOWLER/ised/*')

   else:
      ##  Single stellar population.
      for metallicity in ['m42', 'm62']:
         files.append(root + '/bc2003_hr_%s_salp_ssp.ised' % metallicity)
         
   for fname in files:
     ##  SED ordered by increasing frequency [ergs/s/Hz], Age [GYr], (Monotonically increasing) vs [Hertz]; (Monotonically decreasing) Waves [AA].                                                                       
     seds, ages, vs, ls = read_ised(fname)

     if alliseds:
        metallicity, imf, dummy, xx, xxx = fname.split('bc2003_lr_')[1].split('.ised')[0].split('_')
        tau                              = np.float(xx + '.' + xxx)

     else:
        metallicity = fname.split('_')[3]
        imf         = 'salp'
        tau         = np.NaN

     for ii, sed in enumerate(seds.T):
       if (ages[ii] > 0.01) & (ii % 25 == 0):
         ##  Sample of Ages and metallicities.
         ##  vs * Fv = ls * Fl.  Fig. 13.5 of Cosmological Physics.                                                                                          

         if metallicity == 'm42':
           Ms.append(0.5)  ##  Zo
         
         elif metallicity == 'm62':
           Ms.append(1.0)  ##  Zo
         
         else:
           raise UserWarning('Hard code metallicity for given BC03 mXX: %s' % metallicity)

         Is.append(imf)
         Ts.append(tau)  
         Ages.append(ages[ii])
         spectra.append(sed)

   nmax = len(redshifts) * len(spectra) * len(magnitudes)

   if (ngal is None):
     ngal =  nmax

   elif(ngal > nmax):
     ngal =  nmax

     raise  ValueWarning('Requested LBG number is larger than that available.')

   else:
     pass
   
   runs = [[z, specid, mag] for mag in magnitudes for z in redshifts for specid in np.arange(len(spectra))]
   
   if test:
     runs    =  [[2.3, 0, 22.], [2.3, 1, 22.], [1.7, 2, 22.], [1.7, 3, 22.], [2.2, 4, 22.]]
     printit =  True
     ngal    =  5

   if restframe:
     ngal       =  len(spectra)
     runs       =  [[0.0, specid, 23.5] for specid in np.arange(len(spectra))]
     printit    =  True

     print(runs)

   ##  nmax random selection of galaxies.
   runs = np.array(runs)
   runs = runs[np.random.choice(runs.shape[0], ngal, replace=False)]
   
   ##  flux_factors = 10.0 ** (0.4 * (magnitudes[0] - magnitudes))
   for run in runs:
     redshift  =  np.float(run[0])
     specid    =   np.rint(run[1]).astype(np.int)
     magnitude =  np.float(run[2])

     EBV       =  0.1 * choice(5, 1, replace=True)

     if printit:
       print('\nSolving for redshift: %.3lf, mag: %.3lf, spec. Age: %s, Spec. Metallicity: %s, EBV: %.3lf' % (redshift, magnitude, str(Ages[specid]), Ms[specid], EBV)) 

     ##  Spectra wavelength is very large.  Cut to 2.e4 AA.
     waves   =  ls[ls < 2.e4]
     Fv      =  spectra[specid][ls < 2.e4]
     
     if calzetti:
       Av    =  4.05 * EBV 
       
       waves =  waves.astype(np.double)
       Fv    =  ext_apply(calzetti00(waves, a_v=Av, r_v=4.05, unit='aa'), Fv)

     if madau & (not restframe):
       Fv    =  lephare_madau(waves[::-1], redshift)[::-1] * Fv

     waves   =  waves * u.AA
     Fv      =  Fv * u.erg / u.s / u.cm / u.cm / u.Hz

     Fl      =  Fv.to(u.erg / u.s / u.cm / u.cm / u.AA, equivalencies = u.spectral_density(waves)) 

     ##  Change wavelength sampling to linear with dlambda = 0.1A.
     xx      =  np.arange(waves.value.min(), waves.value.max(), 0.1)

     ##  Note:  reverse to monotically increasing wavelength for resample_flux.
     Fl      =  resample_flux(xx, waves.value[::-1], Fl.value[::-1]) * u.erg / u.s / u.cm / u.cm / u.AA
     waves   =  xx * u.AA

     Fv      =  Fl.to(u.erg / u.s / u.cm / u.cm / u.Hz, equivalencies = u.spectral_density(waves))
     vs      =  waves.to(u.Hz, equivalencies = u.spectral())

     if not restframe:
       ## Redshift spectra for quickspectra input.  'Else', leave unredshifted for redrock template.
       vs, Fv  =   redshift_spectra(vs.value, Fv.value, redshift)

       vs      =   vs * u.Hz
       waves   =   vs.to(u.AA, equivalencies = u.spectral())     ##  Redshifted wavelength.

       Fv      =   Fv * u.erg / u.s / u.cm / u.cm / u.Hz  

     mags      =   get_appmags(vs.value, Fv.value, filters, printit = printit)

     dband     =   get_detband(redshift, lsstfilters, printit=False, restframe=restframe)
     Fv       *=   10. ** (0.4 * (mags[dband] - magnitude))      ##  Normalise to magnitudes[0], e.g. 20.0 in g.

     ##  And check.
     ##  mags  = get_appmags(vs, Fv, filters, printit = True)
       
     ##  F_lambda from F_nu.
     Fl  = Fv.to(u.erg / u.s / u.cm / u.cm / u.AA, equivalencies = u.spectral_density(waves))       ##  [       ergs / s / sq. cm / A]
       
     if restframe:
       wave  = waves.value
       Fl   /= 1.e-17  

       owave = xx

     else:
       Fl   /= 1.e-17                                                                                ##  [1.e-17 ergs / s / sq. cm / A]        
       Fl    = resample_flux(owave, waves.value, Fl.value) * u.erg / u.s / u.cm / u.cm / u.AA

     ##  Magnitude is arbitrary for redrock template. 
     objtypes.append('BC03')
     
     if alliseds:
       subtypes.append('ALL')
     
     else:
       subtypes.append('SSP')

     meta_zs.append(redshift)
     meta_dbands.append(dband)
     meta_mags.append(magnitude)
     meta_ages.append(Ages[specid]) 
     meta_metals.append(Ms[specid])
     meta_imfs.append(Is[specid])
     meta_taus.append(Ts[specid])
     meta_calz.append(EBV)

     flux.append(Fl.value)

     count += 1

     if count >= ngal:
       flux =  tuple(flux)
       flux =  np.column_stack(flux)
       flux =  flux.T

       meta =  Table([objtypes, subtypes, meta_zs, meta_ages, meta_metals, meta_imfs, meta_taus, meta_calz, meta_dbands, meta_mags],\
                      names=('OBJTYPE', 'SUBTYPE', 'REDSHIFT', 'AGE', 'METALLICITY', 'IMF', 'Tau', 'CalzEBV', 'dband', 'AB mag'),\
                       meta={'name': 'BC03 MAKER META'})

       if flux.min() < 0.0:
         raise ValueError('ERROR:  z: %.3lf, spectra: %s' % (redshift, specname))

       return  flux, owave, meta
   

if __name__ == '__main__':
  from  utils  import  latexify


  ##  latexify(fig_width=None, fig_height=None, columns=2, equal=False)

  print('\n\nWelcome to (BC03) LBG maker.')

  ngal              =      15    ##  Over written if rest-frame.  
  test              =   False
  alliseds          =    True
  calz              =    True    ##  Calzetti ext.
  madau             =    True    ##  Madau ext.
  target_type       =   'BC03'
  save              =   False
  restframe         =   False

  redshifts         =   2.0 + np.linspace(0.,   4.0,    20)
  magnitudes        =  20.0 + np.arange(0.,   7.0,   0.1)

  flux, wave, meta  =  BC03_maker(ngal=ngal, restframe=restframe, printit=False, test=test, redshifts=redshifts, magnitudes=magnitudes, alliseds=alliseds, calzetti=calz, madau=madau)

  uwave             =  wave * u.AA 
  vs                =  uwave.to(u.Hz, equivalencies=u.spectral())

  print
  print(wave.shape)
  print(flux.shape)
  print(meta)

  exit(0)

  if restframe & (ngal <= 5):
    filters   = prep_filters(['STEIDEL', 'SUBARU', 'JKC', 'HUBBLE'])

    ##  UGVRI.
    for i, band in enumerate(['U', 'G', 'V', 'R', 'I']):
      if filters[band]['ppkey'][-2] is not 'y':
          pl.fill_between(filters[band]['ls'], filters[band]['Ts'] * 10., alpha=0.3, label='', color='k')   ##  label=filters[band]['ppkey']
    
    for i, band in enumerate(['B']):
      if filters[band]['ppkey'][-2] is not 'y':
          pl.fill_between(filters[band]['ls'], filters[band]['Ts'] * 5., alpha=0.4, label='')

    for i, band in enumerate(['acs_f435w', 'acs_f606w', 'acs_f775w', 'acs_f850lp']):
      band = band.upper()

      if filters[band]['ppkey'][-2] is not 'y':
          pl.fill_between(filters[band]['ls'], filters[band]['Ts'] * 1., alpha=0.4, label='')

    ##  Ages sort. 
    ages   = meta['AGE'].quantity
    ainds  = ages.argsort()

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    count  = 0

    for i, x in enumerate(flux[ainds,:]):
      interim  =        x * u.erg / u.s / u.cm / u.cm / u.AA
      Fv       = interim.to(u.erg / u.s / u.cm / u.cm / u.Hz, equivalencies = u.spectral_density(wave * u.AA))

      if meta['METALLICITY'][ainds][i] == 1.0:
        if count in [1, 2, 3, 4]:
          pl.loglog(wave, wave * interim, label=r'%.2lfGyr' % meta['AGE'][ainds][i], color=colors[count])
  
        count += 1

      else:
        if count in [1, 2, 3, 4]:
          pl.loglog(wave, wave * interim, label=r'', color=colors[count], alpha=0.4)

    pl.xlim(4.e2,  1.1e4)
    pl.ylim(1.e-1, 1.0e4)

    pl.xlabel(r'$\lambda \ [\AA]$')
    pl.ylabel(r'$\lambda \cdot F_\lambda \ [10^{-17} \ \rm{ergs} / s / \rm{cm}^2]$')

    pl.legend(ncol=2, loc=3)

    pl.show(block=True)
    ##  pl.savefig(root + '/gal_maker/plots/BC03.pdf')

    exit(0)

  elif ngal <= 5:
    filters   = prep_filters(['STEIDEL', 'SUBARU', 'JKC', 'HUBBLE'])

    for i, band in enumerate(filters.keys()):
      if filters[band]['ppkey'][-2] is not 'y':
          pl.fill_between(filters[band]['ls'], filters[band]['Ts'] * 0.85, alpha=0.4, label=filters[band]['ppkey'])

    for i, x in enumerate(flux):
      interim =        x * u.erg / u.s / u.cm / u.cm / u.AA
      Fv      = interim.to(u.erg / u.s / u.cm / u.cm / u.Hz, equivalencies = u.spectral_density(wave * u.AA))

      pl.semilogy(wave, vs * Fv, label=str(meta[dband][i]) + ', %.2lfGy, %.2lf' % (meta['AGE'][i], meta['REDSHIFT'][i]))

    pl.xlim(3.6e3, 1.e4)
    pl.ylim(1.e-1, 1.e5)

    pl.xlabel(r'Wavelength [Angstroms]')

    pl.legend(ncol=2, loc=4)
    pl.show(block=True)

  else:
    pass
  
  if save:
    meta.write(os.environ['BEAST'] + '/gal_maker/dat/Tables/galmaker-%s-meta.txt' % target_type, format='ascii', overwrite=True)

    simdir = os.path.join(os.environ['SCRATCH'], 'desi', 'simspec')
  
    ## os.makedirs(simdir, exist_ok=True)

    infile = os.path.join(simdir, '%s-input-spectra.fits' % target_type)

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
