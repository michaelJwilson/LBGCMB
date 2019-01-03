import numpy             as     np
from   scipy.interpolate import interp1d


def add_desinoise(ls, Fl):
    """
    Place check that dlambda of ls is the same                                                                                                                                                                                           
    as in the noise curve.
    """

    noise       = np.loadtxt("../quickspectra/nosource-exp1800-singlearm.txt")  ## [1e-17 ergs/s/cm2/A]
    noise[:,1] *= 1.e-17
    
    interp      = interp1d(noise[:,0], noise[:,1], bounds_error=False, fill_value=0.0) 

    Fl         += interp(ls)

    return ls, Fl
