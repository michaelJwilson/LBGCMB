from    __future__           import  absolute_import, division, print_function

import  time 
import  os 
import  sys
import  time
import  argparse
import  collections

import  desitarget
import  desispec.io
import  desimodel.io
import  desisim.io
import  desisim.obs
import  desisim.simexp

import  numpy                as      np
import  astropy.table
import  astropy.units        as      u
import  astropy.io.fits      as      pyfits

from    desiutil.log         import  get_logger
from    desispec.spectra     import  Spectra, spectra_dtype
from    desispec.resolution  import  Resolution


def sim_spectra(wave, flux, program, spectra_filename, obsconditions = None, expid = 0, seed = 0, survey = 'desi'):
    '''
    Simulate spectra from input observer wavelength and (redshifted) flux and writes a .FITS file in the spectra 
    format that can be used as input to the redshift fitter.

    Args:
        Wave : 1D np.array of wavelength in Angstrom (in vacuum) in observer frame (i.e. redshifted)
        Flux : 1D or 2D np.array. 1D array must have same size as wave, 
               2D array must have shape[1] = wave.size for multiple input. 
               
               Note:  Flux has to be in units of 1e-17 [ergs/s/cm2/A].

        spectra_filename:  Path to output FITS file in the Spectra format
    
    Optional:
        obsconditions:  Dictionary of observation conditions: {SEEING, EXPTIME, AIRMASS, MOONFRAC, MOONALT, MOONSEP}
        expid:          This expid number will be saved in the spectra fibermap
        seed:           Random seed       
    '''

    log = get_logger()

    if len(flux.shape) == 1:
      flux  = flux.reshape((1, flux.size))

    nspec   = flux.shape[0]
    
    log.info("Starting simulation of {} spectra".format(nspec))
    
    tileid  = 0
    
    telera  = 0
    teledec = 0    
    
    night   = desisim.obs.get_night(utc = time.gmtime())
               
    frame_fibermap                = desispec.io.fibermap.empty_fibermap(nspec)    
    frame_fibermap.meta["FLAVOR"] = "custom"
    frame_fibermap.meta["NIGHT"]  = night
    frame_fibermap.meta["EXPID"]  = expid
    
    ##  Add DESI_TARGET and TARGETID
    tm  = desitarget.desi_mask            ## Contains 'templates' for STD_FSTAR.
    
    for spec in range(nspec):
      frame_fibermap['DESI_TARGET'][spec] = tm.STD_FSTAR
      frame_fibermap['TARGETID'][spec]    = spec
        
    ## Spectra fibermap has two extra fields: night and expid.
    spectra_fibermap = np.zeros(shape=(nspec,), dtype=spectra_dtype())
    
    for s in range(nspec):
      for tp in frame_fibermap.dtype.fields:
        spectra_fibermap[s][tp]  = frame_fibermap[s][tp]
    
    spectra_fibermap[:]['EXPID'] = expid  ## Needed by spectra.
    spectra_fibermap[:]['NIGHT'] = night  ## Needed by spectra.

    program = program.lower()
    
    if obsconditions is None:
      if program in ['dark', 'lrg', 'qso']:
        """
        E.g.
          reference_conditions['DARK']['SEEING']     =  1.1
          reference_conditions['DARK']['EXPTIME']    = 1000
          reference_conditions['DARK']['AIRMASS']    =  1.0
          reference_conditions['DARK']['MOONFRAC']   =  0.0
          reference_conditions['DARK']['MOONALT']    =  -60
          reference_conditions['DARK']['MOONSEP']    =  180
        """
        obsconditions  = desisim.simexp.reference_conditions['DARK']
    
      elif program in ['elg', 'gray', 'grey']:
        obsconditions  = desisim.simexp.reference_conditions['GRAY']
        
      elif program in ['mws', 'bgs', 'bright']:
        obsconditions  = desisim.simexp.reference_conditions['BRIGHT']

      else:
        raise  ValueError('Unknown program {}'.format(program))
    
    elif isinstance(obsconditions, str):
      try:
        obsconditions  = desisim.simexp.reference_conditions[obsconditions.upper()]
        
      except KeyError:
        raise ValueError('Input observation conditions {} are not in {}'.format(obsconditions.upper(), list(desisim.simexp.reference_conditions.keys())))
    
    if  survey == 'pfs':
        wavemin   =   3796.  ## [A] 
        wavemax   =  12605.  ## [A]

    elif  survey == 'desi':    
        wavemin   = desimodel.io.load_throughput('b').wavemin
        wavemax   = desimodel.io.load_throughput('z').wavemax

    elif  survey == 'beast':
        wavemin   =   3796.  ## [A]
        wavemax   =  12605.  ## [A]

    else:
        raise ValueError("\n\nSurvey %s is not available." % survey)
    
    log.info("Setting wave limits to survey: {} ... {} to {} [A]".format(survey, wavemin, wavemax))
    
    ii           = (wavemin <= wave) & (wave <= wavemax)

    flux_unit    = 1e-17 * u.erg / (u.Angstrom * u.s * u.cm ** 2)
    
    wave         =   wave[ii] * u.Angstrom  ## Dimensionful quantities; only between wavelength limits.  
    flux         = flux[:,ii] * flux_unit
    
    sim          = desisim.simexp.simulate_spectra(wave, flux, fibermap = frame_fibermap, obsconditions = obsconditions, survey=survey)
    
    ## Add random noise. 
    random_state = np.random.RandomState(seed)

    sim.generate_random_noise(random_state)
    
    specdata     = None
    scale        = 1.e17

    resolution   = {}
    
    ## Methods for sim object. 
    for camera in sim.instrument.cameras:
        R                       = Resolution(camera.get_output_resolution_matrix())
        resolution[camera.name] = np.tile(R.to_fits_array(), [nspec, 1, 1])
        
    for table in sim.camera_output:        
        wave =  table['wavelength'].astype(float)
        flux = (table['observed_flux'] + table['random_noise_electrons'] * table['flux_calibration']).T.astype(float)
        ivar =  table['flux_inverse_variance'].T.astype(float)
        
        band =  table.meta['name'].strip()[0]
        
        flux =  flux * scale
        ivar =  ivar / scale**2
        mask =  np.zeros(flux.shape).astype(int)
        
        ## Create spectra object for redrock. 
        spec = Spectra([band], {band : wave}, {band : flux}, {band : ivar}, 
                       resolution_data = {band: resolution[band]}, 
                       mask            = {band: mask}, 
                       fibermap        =  spectra_fibermap, 
                       meta            =  None,
                       single          =  True)
        
        if specdata is None:
            specdata = spec

        else:
            specdata.update(spec)
    
    log.info("Writing to: %s" % spectra_filename)

    desispec.io.write_spectra(spectra_filename, specdata)        

    log.info('Successfully created %s.' % spectra_filename)

def get_input(input, repeat=1):
    ## Get input file.                                                                                                                                     
    try:
      hdulist         = pyfits.open(input)

      input_wave      = hdulist["WAVELENGTH"].data
      input_flux      = np.file(hdulist["FLUX"].data, (repeat, 1))

      hdulist.close()

    except(OSError, IOError):
      try:          
          tmp         = np.loadtxt(input).T
          
          input_wave  = tmp[0]                         ## First column assumed to be wavelength.                                                      
          input_flux  = np.tile(tmp[1:], (repeat, 1))  ## Remaining columns assumed to be input spectra.                                            

          log.info("Input flux shape (after repeat) = {}".format(input_flux.shape))
          
      except(ValueError, TypeError):
          log.error("Could not read ASCII file; Needs at least two columns separated by ' ',")
          log.error("The first is wavelength [A in vacua], the others for flux [1e-17 ergs/s/A/cm2], one column per spectrum.")
          
          sys.exit(12)
          
      return  input_wave, input_flux


if __name__ == "__main__":
    print("\n\nWelcome to quickspectra.\n\n")

    log                         =  get_logger()

    seed                        =  0

    ##  survey                  = 'desi' 
    ##  survey                  = 'pfs'
    survey                      = 'beast'
    

    program                     = 'DARK'

    start_time                  =  time.time()
    sample_rate                 =            1

    ##  BC03 production run.                                                                                                                               
    ##  redshifts    =  3.5 + np.linspace(0., 2.50, 25)[::sample_rate]
    ##  magnitudes   = 20.0 + np.linspace(0., 6.00, 25)[::sample_rate]

    ##  exposures    =   60. * 15.* np.arange(1, 22, 1)[::sample_rate]

    ##  Shapley composite production run;
    redshifts   =   1.5 + np.linspace(0., 1.50, 15)[::sample_rate]                ##  3.5 + np.linspace(0., 2.50, 25)
    magnitudes  =  20.0 +   np.arange(0., 7.00,  1)                               ##  Magnitudes is implicitly looped over.  QuickSpectra generates 
                                                                                  ##  noise curves for each column in input file.
    exposures    =  60. * 15. * np.arange(1, 30, 3)[::sample_rate]                ##  Seconds. 

    print('\n\nSolving for redshifts:\n'  + ''.join('%.3lf, ' % x for x in redshifts))
    print('\n\nSolving for magnitudes:\n' + ''.join('%.3lf, ' % x for x in magnitudes))
    print('\n\nSolving for exposures:\n'  + ''.join('%.3lf, ' % x for x in exposures))

    obsconditions               = desisim.simexp.reference_conditions[program]

    nruns                       = len(exposures) * len(redshifts)

    for ii, exptime in enumerate(exposures[:1]):          
      ##  Successfully created /Users/M.J.Wilson/work/LBGCMB//quickspectra/dat/quickspectra/desi/Shapley/Q3/spec-Shapley-Q3_199-z2.375_exp12600.fits.  
      obsconditions['EXPTIME']  = exptime  

      for jj, redshift in enumerate(redshifts[:1]):
        ##  input  = os.environ['LBGCMB'] + '/quickspectra/dat/in_quickspectra/BC03/spec-BC03-z%.1lf.dat'           %   redshift
        ##  input  = os.environ['LBGCMB'] + '/quickspectra/dat/in_quickspectra/Shapley/spec-Shapley-all-z%.3lf.dat' %   redshift
        input      = os.environ['LBGCMB'] + '/quickspectra/dat/in_quickspectra/baseline/spec-all-z%.3lf.dat'        %   redshift

        ##  output = os.environ['LBGCMB'] + '/quickspectra/dat/quickspectra/%s/BC03/spec-BC03-z%.1lf_exp%.0lf.fits'               %  (survey, redshift, exptime)
        ##  output = os.environ['LBGCMB'] + '/quickspectra/dat/quickspectra/%s/Shapley/all/spec-Shapley-all-z%.3lf_exp%.0lf.fits' %  (survey, redshift, exptime)
        output     = os.environ['LBGCMB'] + '/quickspectra/dat/quickspectra/%s/baseline/spec-all-z%.3lf_exp%.0lf.fits'            %  (survey, redshift, exptime)

        input_wave, input_flux  = get_input(input)
        
        sim_spectra(input_wave, input_flux, input, obsconditions = obsconditions, spectra_filename = output, seed = seed, survey = survey)
        
        print('\n\n\nTime: %.2lf s.  Percentage complete: %.2lf;  Redshift: %.2lf, Exposure: %.2lf \n\n\n' % (time.time() - start_time,\
                                                                                                              100. * (jj + ii * len(redshifts)) / nruns,\
                                                                                                              redshift, exptime))
    
    print("\n\nDone.\n\n")
