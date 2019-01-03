import  numpy              as      np
import  astropy.units      as      u 
import  astropy.constants  as      const
import  itertools          as      it

from    prep_filters       import  prep_filters
from    redshift_spectra   import  redshift_spectra
from    app_mags           import  get_appmags
from    read_ised          import  read_ised
from    scipy.interpolate  import  interp1d 
from    rrtemplate_io      import  create_template


if __name__ == "__main__":
  print('\n\nCreate input curves for quickspectra from\nBC03 templates.\n\n')

  template              =  True
    
  ##  Age [Myr]; vs [Hz]; SEDS [ergs/s/Hz]; ls [A]; Ll [ergs/s/angstrom].                                                                                   
  ##  Ordered by wavelength.
  ages, vs, Fv, ls, Ll  =  read_ised('GALAXEV/models/Padova1994/salpeter/bc2003_hr_m72_salp_ssp.ised', AgeMyr = 25., printit = False)
  
  ##  Wavelengths for output.
  waves                 =  np.arange(1., 1.5e4, 1.)

  ## Prep. filters.                                                                                                                                       
  filters               =  prep_filters()

  if template:
    mags                =  get_appmags(vs, Fv, filters, printit = False)
    Ll                 *=  10. ** (0.4 * mags['i'] - 0.4 * 22.5)

    owaves              =  np.arange(1., 1.e4, 0.1)
    interp              =  interp1d(ls.value, Ll.value, kind='linear', bounds_error=False, fill_value=0.0, assume_sorted=False)

    ##  np.savetxt('../quickspectra/spectra/BC03/restframe/spec-BC03-z%.1lf.dat' % 0.0, np.c_[waves, interp(waves)],\
    ##             fmt='%.6le', header='   WAVELENGTH        FLUX\n--------------------------')

    ##  Twinned with unredshifted spectra for redrock template. .fits output.                                                                              
    ##  create_template(owaves, interp(owaves) / 1.e-17, 'LBG', printit = False)
    create_template(ls.value, Ll.value / 1.e-17, 'LBG', printit = False)   

    exit()

  ##  Normalise to -250th mag. as baseline; problem with underflow.                                                                                         
  mags        =  get_appmags(vs, Fv, filters, printit = False)
  
  Fv         *=  10. ** (0.4 * mags['i'] - 0.4 * -250.0)

  ##  Set up with baseline vs and Fv.  Now renormalise and redshift spectra:                                                                             
  redshifts   =   2.0 + np.arange(0., 2.60, 0.10)
  magnitudes  =  20.0 + np.arange(0., 6.00, 0.25)

  result      =  [waves]

  for redshift in redshifts:
    ovs, oFv  =  redshift_spectra(vs, Fv, redshift)  
     
    ols       =  ovs.to(u.AA, equivalencies = u.spectral())

    ##  Order by increasing wavelength.
    ols       =  ols[::-1]

    ##  Magnitudes of redshifted Fv.
    omags     =  get_appmags(ovs, oFv, filters, printit = False)
        
    for magnitude in magnitudes:  
      print("Calculating z = %.3lf, m = %.3lf." % (redshift, magnitude))

      ## The quickspectra input file is an ASCII file with wavelength in [A] as                                                                            
      ## the first column (in vacuum, REDSHIFTED); the other columns are treated                                                                           
      ## as spectral flux densities in units of 1e-17 [ergs/s/cm2/A].                                                                                      
      mFv      =  oFv * 10. ** (0.4 * omags['i'] - 0.4 * magnitude)
      mFv     /=  1.e-17

      ## mmags =  get_appmags(ovs, mFv, filters, printit = False)

      Ll       = (mFv * ovs**2.) / const.c.value               ## [ergs/s/meter].                                                                           
      Ll      *=  1e-10                                        ## [ergs/s/angstrom]. 

      ##  Order by increasing wavelength;  Previously dropped.
      Ll       =  Ll[::-1]

      interp   =  interp1d(ols.value, Ll.value, kind='linear', bounds_error=False, fill_value=0.0, assume_sorted=False)

      ## Input to quick spec. 
      result.append(interp(waves))
         
    np.savetxt('../quickspectra/dat/in_quickspectra/BC03/spec-BC03-z%.1lf.dat' % redshift, np.array(result).T,\
               fmt="%.12le", header='   WAVELENGTH        FLUX\n--------------------------')
  
  print("\n\nDone.\n\n")

