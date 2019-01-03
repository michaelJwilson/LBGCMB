import  os
import  glob
import  numpy               as      np
import  pylab               as      pl
import  astropy.constants   as      const
import  astropy.units       as      u

from    scipy.interpolate   import  interp1d
from    prep_filters        import  prep_filters
from    read_ised           import  read_ised
from    rrtemplate_io       import  read_template
from    utils               import  luminosity_distance 


def load_star(type=None):
  """                                                                                                                                                     
  Load ls, Fl of stellar type.
  """

  if type is None:
    type    = 'b5iii'

  ## 
  fname     = os.environ['LBGCMB'] + '/dropouts/stellar_spectra/PICKLES/%s.sed' % str(type)

  spectra   = np.loadtxt(fname)

  ls        = spectra[:,0]
  Ll        = spectra[:,1]

  return  ls, Ll

def get_appmags(vs, Fv, filters, printit=False): 
  """
  Given a list of filters, rest wavelengths ls and F_\lambda at those wavelengths, 
  return the apparent magnitude in each filter band.
  """
  
  import  collections

  iFv            = interp1d(vs, Fv, bounds_error=False, fill_value=0.0)                      ## Interpolate L_lambda.
  
  mags           = collections.OrderedDict()

  ## Implement eqn. (7) of ./GALAXEV/doc/bc03.pdf 
  for i, band in enumerate(filters.keys()):                                         
    filter       = filters[band]

    ## Filters are defined in wavelength; calculate normalisation by eqn. (8) denominator.
    norm         = np.trapz(filter['Ts'] / filter['vs'], filter['vs'])

    result       = np.trapz(iFv(filter['vs']) * filter['Ts'] / filter['vs'], filter['vs'])
    result      /= norm

    mm           = -2.5*np.log10(result) - 48.60                                             ## AB bandpass magnitude.

    mags[band]   = mm

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


if __name__ == "__main__":
  print("\n\nWelcome to app_mags.\n\n")

  filters                    = prep_filters()

  types                      = glob.glob(os.environ['LBGCMB'] + 'dropouts/stellar_spectra/PICKLES/*')
  types                      = [type.split('/')[-1].split('.sed')[0] for type in types] 
  
  LBGCOLORS                  = ['u-g', 'g-r', 'r-i', 'i-z']

  for type in types:
    ls, Fl                   = load_star(type)                         

    print("\n\nCalculating magnitudes for type: %s.\n\n" % type)

    vs                       = (1e10 / ls) * const.c.to('m/s').value                            
    Fv                       = ls * Fl / vs
    
    fname                    = '../colors/colors_%s.json' % str(type)

    mags                     = get_appmags(vs, Fv, filters, printit =  True)
    colors                   =  get_colors(mags,  LBGCOLORS,  fname = fname)
  
  print("\n\nDone.\n\n")
