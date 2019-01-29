import  numpy               as      np
import  pylab               as      pl
import  astropy.constants   as      const
import  astropy.units       as      u

from    scipy.interpolate   import  interp1d
from    prep_filters        import  prep_filters
from    read_ised           import  read_ised
from    rrtemplate_io       import  read_template
from    utils               import  luminosity_distance 
from    composite           import  load_composite
from    magABsource         import  magAB
from    redshift_spectra    import  redshift_spectra
from    uvlum               import  get_uvlum
from    add_desinoise       import  add_desinoise


def load_sed(fname="../quickspectra/spectra/spec-ABmag22.0.dat", flux_unit=1.e-17):
  '''                                                                                                                                                  
  Load ls, Fl of AB magnitude source of constant Fv;                                                                                                  
  Default is 22 mag, corresponding to Fv of 5.75 micro Jy.                                                                                               
  '''

  spectra   = np.loadtxt(fname)

  ls        = spectra[:,0]
  Ll        = spectra[:,1]

  return ls, Ll * flux_unit


def get_appmags(vs, Fv, filters, printit=False): 
  '''
  Given a list of filters, rest wavelengths ls and F_\lambda at those wavelengths, 
  return the apparent magnitude in each filter band.
  '''
  
  import  collections

  iFv            = interp1d(vs, Fv, bounds_error=False, fill_value=0.0)          ## Interpolate L_lambda.
  
  mags           = collections.OrderedDict()

  ##  Check vs is monotonically increasing.                                                                                                                
  assert  np.all(np.diff(vs) >= 0.)

  ## Implement eqn. (7) of ./GALAXEV/doc/bc03.pdf 
  for i, band in enumerate(filters.keys()):                                         
    filter       = filters[band]
    
    ##  Assume filter vs are derived from ls and are monotonically decreasing,                                                                             
    ##  in which case they need reversed.                                                                                                                    
    assert np.all(np.diff(filter['vs']) <= 0.)

    ## Filters are defined in wavelength; calculate normalisation by eqn. (8) denominator.
    norm         = np.trapz(filter['Ts'] / filter['vs'], filter['vs'])

    result       = np.trapz(iFv(filter['vs']) * filter['Ts'] / filter['vs'], filter['vs'])
    result      /= norm

    if result == 0.0:
      mags[band] = -1. * np.NaN

    else:
      mags[band] = -2.5 * np.log10(result) - 48.60                                           ##  AB bandpass magnitude.           

    if printit:
      print("%s \t %.6f" % (band, mags[band]))

  return  mags

def get_colors(mags, get_colors=None, fname = None):
  """
  Given magnitudes as input, return either
  an OrderedDict of all possible colors, or 
  the colors requested. 
  """
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

def plot_Fv(filters, ls, Fv):
  pl.plot(ls, Fv)

  for i, band in enumerate(filters.keys()):
    pl.plot(filters[band]['ls'], filters[band]['Ts'] * Fv[-1], label = filters[band]['ppkey'])

  pl.xlim(3.e3, 1.e4)
  pl.ylim(-Fv[-1] / 10., 1.4 * Fv[-1])

  pl.xlabel(r'Wavelength [$\AA$]')
  pl.ylabel(r'$F_{\nu}$')

  pl.savefig("plots/mag22AB_LyLimBreak.pdf")


if __name__ == "__main__":
  print("\n\nWelcome to app_mags.\n\n")

  filters                    = prep_filters()

  ## BStar                   = read_template(type='star-B', printit=False)

  ## ls, Fl                  = load_sed("./dat/elg-avg-line-ratio.dat")                          ## ls [AA], Ll [ergs/s/A].
  ls, Fl                     = load_sed("./dat/spec-lrg-z0.7-zmag20.00.dat")                     ## ls [AA], Ll [ergs/s/A].   

  print ls
  print Fl

  ## ls, Fv                  = load_composite(mJy=False, getcontinuum=False, wext=1.2e4)
  ## ages, vs, Lv, ls, Ll    = read_ised('GALAXEV/models/Padova1994/salpeter/bc2003_hr_m72_salp_ssp.ised', AgeMyr=25., printit = False)

  LBGCOLORS                  = ['u-g', 'g-r', 'r-i', 'i-z']

  ## lbg_redshifts           = np.arange(1.0,  7.5, 0.125)
  lrg_redshifts              = np.arange(0.00, 1.3, 0.050)
  
  for redshift in lrg_redshifts:  
    print("\n\nCalculating magnitudes at z = %.3lf.\n\n" % redshift)

    dlambda                  = 1.0 

    ## ls                    = np.arange(0., 1.e4, dlambda)                                      ## Restframe [A].
    ## vs                    = (1e10 / ls) * const.c.to('m/s').value                             ## ls [AA]; vs [Hz]; c [m/s.]
    
    ## ls                    = (1e10 / vs) * const.c.to('m/s').value                             ## Redshifted wavelengths.  
    ## Fl                    = vs * Fv / ls

    ## ls, Fl                = add_desinoise(ls, Fl)

    vs                       = (1e10 / ls) * const.c.to('m/s').value                              
    Fv                       =  ls * Fl / vs

    ## Luv                   = get_uvlum(ls, Ll) 
    
    ## vs, Fv                = magAB(vs, 25.0, redshift = redshift, line = 912., depth=26.5)     ## Redshifts frequencies.                                 
    
    vs, Fv                   = redshift_spectra(vs, Fv, redshift)                                ## redshifted spectra, observed frequencies. 

    fname                    = 'colors/LRG/colors_lrg_z%.2lf.json' % redshift

    mags                     = get_appmags(vs, Fv, filters, printit =   True)
    colors                   =  get_colors(mags,  LBGCOLORS,  fname =  fname)

    ## if redshift == 4.0:
    ##   plot_Fv(filters, ls, Fv)
  
  print("\n\nDone.\n\n")
