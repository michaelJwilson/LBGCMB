import  os
import  numpy             as      np

from    specs             import  samplestats
from    cosmo             import  cosmo
from    scipy.interpolate import  interp1d


def galaxy_frac(ms, dropband='g'):
    if not dropband in ['g', 'r']:
        raise  ValueError('\n\nGalaxy fraction not available for dropband: %s' % dropband)

    ##  Bottom panel of Fig. 7 top of https://arxiv.org/pdf/1704.06004.pdf;                                                                                 
    ##  Galaxy fraction of GOLDRUSH g-dropouts.                                                                                                             
    fname = os.environ['LBGCMB'] + '/dropouts/goldrush/cats/galaxy_fraction.dat'
    data  = np.loadtxt(fname)

    if dropband == 'g':
      frac  = interp1d(data[:,0], data[:,1], kind='linear', bounds_error=False, fill_value=(0.0, 1.0))
                       
    else:
      frac  = interp1d(data[:,0], data[:,2], kind='linear', bounds_error=False, fill_value=(0.0, 1.0))
             
    ##  Conversion from Ms to ms, e.g. eqn (13) of https://arxiv.org/pdf/1704.06004.pdf                                                              
    specs = samplestats()

    zeff  = specs[dropband]['z']

    Ms    = ms + 2.5 * np.log10(1. + zeff) - 5. * np.log10(1.e6 * cosmo.luminosity_distance([zeff]).value / 10.)

    return  frac(Ms)
