import  os
import  numpy                   as      np
import  pylab                   as      pl
import  astropy.constants       as      const
import  astropy.units           as      u

from    scipy.interpolate       import  interp1d
from    prep_filters            import  prep_filters
from    desispec.interpolation  import  resample_flux
from    magABsource             import  magAB
from    tophat_filters          import  tophat_filters


def get_appmags(vs, Fv, filters, printit=False): 
  '''
  Given a list of filters, rest wavelengths ls and F_\lambda at those wavelengths, 
  return the apparent magnitude in each filter band.
  '''
  
  import  collections


  mags       =  collections.OrderedDict()

  ##  Implement eqn. (7) of ./GALAXEV/doc/bc03.pdf 
  for i, band in enumerate(filters.keys()):                                         
    filter   =  filters[band]

    ##  Note:  Filters extend in wavelength much further than their transmission. 
    ##  assert vs.max() >= filter['vs'].max()
    ##  assert vs.min() <= filter['vs'].min()
    
    pFv      =  resample_flux(filter['vs'][::-1], vs[::-1], Fv[::-1], ivar=None, extrapolate=False)
    pFv      =  pFv[::-1]

    ##  Filters are defined in wavelength;  Calculate normalisation by eqn. (8) denominator.
    norm     =  np.trapz(filter['Ts'], filter['vs'])        ##  / filter['vs']
    result   =  np.trapz(pFv * filter['Ts'], filter['vs'])  ##  / filter['vs']

    result  /=  norm

    if result   == 0.0:
      mags[band] =  99.

    else:
      mags[band] = -2.5 * np.log10(result) - 48.60                                           ##  AB bandpass magnitude.

    if printit:
      print("%s \t %.6le \t %.6lf" % (band, result, mags[band]))
    
  return  mags

def get_colors(mags, get_colors=None, fname = None):
  '''
  Given magnitudes as input, return either
  an OrderedDict of all possible colors, or 
  the colors requested. 
  '''

  import  json
  import  itertools
  import  collections


  colors = collections.OrderedDict()

  if get_colors is None:
    ## Get all possible colors from the magnitudes provided. 
    get_colors = list(itertools.combinations(mags, 2))

  for scolor in get_colors:
    color          = mags[scolor[0]] - mags[scolor[2]] 
    colors[scolor] = color

  if fname is not None:
    print("\n\nWriting colors to %s." % fname)

    json.dump(colors, open(fname, 'w'))    

  return  colors


if __name__ == "__main__":
  from astropy.io import fits


  print("\n\nWelcome to app_mags.\n\n")

  spectra      =  fits.open(os.environ['CSCRATCH'] + '/desi/simspec/safe/bgs-input-spectra.fits')

  ##  filters  =  prep_filters(['LSST'])
  filters      =  tophat_filters()

  wave         =  spectra[0].data * u.AA                                              ##  Angstroms.
  flux         =  spectra[1].data[1,:] * 1.e-17 * u.erg / u.s / u.cm / u.cm /u.AA     ##  1e-17 erg / s / cm2 / Angstrom.

  vs           =  wave.to(u.Hz, equivalencies=u.spectral()).value
  Fv           =  flux.to(u.erg / u.s / u.cm / u.cm / u.Hz, equivalencies=u.spectral_density(wave)).value

  mags         =  get_appmags(vs, Fv, filters, printit=True)

  pl.loglog(vs, Fv)

  print(vs)

  for band in mags:
      ff     =  10. ** -((mags[band]  + 48.60) / 2.5)
      pl.loglog(filters[band]['Cs'], ff, 'c^', markersize=4)

  for band in filters:
    pl.plot(filters[band]['vs'], filters[band]['Ts'])

  pl.ylim(1.e-27, 2.e-25)
  pl.show()

  print("\n\nDone.\n\n")
