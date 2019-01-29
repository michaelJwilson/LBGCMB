import  os
import  numpy                  as      np
import  pylab                  as      pl
import  astropy.constants      as      const
import  astropy.units          as      u

from    scipy.interpolate      import  interp1d
from    prep_filters           import  prep_filters
from    utils                  import  luminosity_distance 
from    magABsource            import  magAB
from    redshift_spectra       import  redshift_spectra
from    desispec.interpolation import  resample_flux


def load_sed(fname, flux_unit=1.e-17 * u.erg / u.s / u.cm / u.cm / u.AA):
  '''                                                                                                                                                     
  Load ls, Fl of AB magnitude source of constant Fv;                                                                                                  
  Default is 22 mag, corresponding to Fv of 5.75 micro Jy.                                                                                               
  '''

  spectra   = np.loadtxt(fname)

  ls        = spectra[:,0] * u.AA
  Fl        = spectra[:,1] * flux_unit

  return  ls, Fl

def get_appmags(vs, Fv, filters, printit=False): 
  '''
  Given a list of filters, rest wavelengths ls and F_\lambda at those wavelengths, 
  return the apparent magnitude in each filter band.
  '''
  
  import  collections
  
  mags           = collections.OrderedDict()

  ## Implement eqn. (7) of ./GALAXEV/doc/bc03.pdf 
  for i, band in enumerate(filters.keys()):                                         
    filter       = filters[band]

    pFv          = resample_flux(filter['vs'][::-1], vs, Fv, ivar=None, extrapolate=False)
    pFv          = pFv[::-1]

    ## Filters are defined in wavelength; calculate normalisation by eqn. (8) denominator.
    norm         = np.trapz(filter['Ts'] / filter['vs'], filter['vs'])

    result       = np.trapz(pFv * filter['Ts'] / filter['vs'], filter['vs'])
    result      /= norm

    if result == 0.0:
      mags[band] = 99.

    else:
      mags[band] = -2.5 * np.log10(result) - 48.60              ##  AB bandpass magnitude. 

    if printit:
      print("%s \t %.6f" % (band, mags[band]))

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
    json.dump(colors, open(fname, 'w'))    

  return  colors

def plot_Fv(filters, ls, Fl, fname='spectra'):
  pl.semilogy(ls, Fl)

  '''
  for i, band in enumerate(filters.keys()):
    pl.plot(filters[band]['ls'], filters[band]['Ts'] * Fv[-1], label = filters[band]['ppkey'])
  '''

  root     =  os.environ['BEAST'] + '/desihub/desimodel/data/spectra/'
  spectra  = ['spec-elg-z0.818.dat', 'spec-elg-z1.189.dat', 'sn-spec-qso-z1.5-rmag23.00.dat', 'spec-lrg-z0.9-zmag20.57.dat']
  spectra  = ['spec-lrg-z0.9-zmag20.57.dat', 'spec-elg-o2flux-8e-17.dat', 'spec-elg-z0.818.dat', 'spec-qso-z1.5-rmag22.81.dat']

  for spec in spectra:
    dat = np.loadtxt(root + spec)
    pl.plot(dat[:,0], dat[:,1], label=spec)

  pl.xlim(3.e3, 1.e4)
  pl.ylim( 0.1,  2.0)

  pl.xlabel(r'Wavelength [$\AA$]')
  pl.ylabel(r'$F_{\lambda}$')

  pl.legend()

  pl.savefig(os.environ['BEAST'] + "/gal_maker/plots/%s.pdf" % fname)


if __name__ == "__main__":
  print("\n\nWelcome to app_mags.\n\n")

  filters                    = prep_filters()
  
  ##  fname                  = os.environ['BEAST'] + '/desihub/desimodel/data/spectra/sn-spec-qso-z2.4-rmag23.00.dat'
  ##  fname                  = os.environ['BEAST'] + '/desihub/desimodel/data/spectra/spec-elg-o2flux-8e-17.dat' 
  fname                      = os.environ['BEAST'] + '/desihub/desimodel/data/spectra/spec-lya.dat'

  LBGCOLORS                  = ['u-g', 'g-r', 'r-i', 'i-z']

  redshifts                  = np.array([1.5])  ##  np.arange(1.0, 1.1, 0.050)
  
  for redshift in redshifts:  
    print("\n\nCalculating magnitudes at z = %.3lf.\n\n" % redshift)

    ls, Fl                   =  load_sed(fname)                                            

    vs                       =  ls.to(u.Hz, equivalencies = u.spectral())
    Fv                       =  Fl.to(u.erg / u.s / u.cm / u.cm / u.Hz, equivalencies = u.spectral_density(ls))
    
    vs, Fv                   =  redshift_spectra(vs.value, Fv.value, redshift)                    ## redshifted spectra, observed frequencies. 
    
    vs                      *=  u.Hz
    ls                       =  vs.to(u.AA, equivalencies = u.spectral())                         ## redshifted wavelengths.  
     
    mags                     =  get_appmags(vs, Fv, filters, printit =  True)
    colors                   =   get_colors(mags,  LBGCOLORS,  fname =  None)

    Fv                       =  Fv * u.erg / u.s / u.cm / u.cm / u.Hz 
    Fv                      *=  10. ** (0.4 * (mags['z'] - 22.81))

    Fl                       =  Fv.to(u.erg / u.s / u.cm / u.cm / u.AA, equivalencies = u.spectral_density(ls))
    Fl                      /=  1.e-17

    mags                     =  get_appmags(vs.value, Fv.value, filters, printit =  True)

    pl.clf()

    plot_Fv(filters, ls.value, Fl.value)
  
  print("\n\nDone.\n\n")
