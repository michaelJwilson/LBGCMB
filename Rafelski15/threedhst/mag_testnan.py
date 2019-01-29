import numpy as np
import astropy.constants as const

from app_mags    import get_appmags
from magABsource import magAB
from prep_filters import prep_filters


magnitude = 26.5
dlambda   =  0.1

ls        = np.arange(3.e3, 1.e4, dlambda)           ## Restframe [A].                                                                                                                                      
vs        = (1.e10 / ls) * const.c.to('m/s').value   ## ls [AA]; vs [Hz]; c [m/s.]                                                                                                                          

vs, Fv    = magAB(vs, magnitude)                     ## Flat Fv        

filters   = prep_filters(['LSST', 'VIDEO', 'STEIDEL', 'SUBARU', 'JKC'], normed=False)
mags      = get_appmags(vs, Fv, filters, printit=True)
