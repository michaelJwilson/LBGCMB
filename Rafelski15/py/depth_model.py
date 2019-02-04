import numpy         as np
import pylab         as pl
import astropy.units as u

from   app_mags      import  get_appmags
from   collections   import  OrderedDict
from   prep_filters  import  prep_filters
from   magABsource   import  magAB


##  Note:  underestimates the noise floor.
##         open question why.  

depths      = OrderedDict()

depths['u'] = 25.8
depths['g'] = 25.4
depths['r'] = 24.9
depths['i'] = 25.0
depths['z'] = 24.7
depths['y'] = 24.4
depths['J'] = 24.0
depths['H'] = 24.0
depths['K'] = 24.0
depths['R'] = 25.5
depths['U'] = 25.5

depths['ACS_F435W']  = 27.0
depths['ACS_F606W']  = 27.0
depths['ACS_F775W']  = 27.0
depths['ACS_F850LP'] = 27.0

band = 'u'
mags = []

ls   = np.arange(3.e3, 1.e4, 0.1) * u.AA
vs   = ls.to(u.Hz, equivalencies=u.spectral()).value

filters = prep_filters(['LSST', 'VIDEO', 'STEIDEL', 'SUBARU', 'JKC', 'HUBBLE'], normed=False)

for x in np.arange(0, 500, 1):
  vs, nFv    =  magAB(vs, mag = depths[band])

  ##  eqn. (10.2) of Chromey, Introduction to Observational Astronomy.                                                                                                                                                             
  planck     =  6.62606885e-27  ## ergs * s.                                                                                                                                                                                       
  lightspeed =  2.9979e10       ## cm/s.                                                                                                                                                                                           

  phi        =  np.rint(vs * nFv / planck / lightspeed).astype(np.int)

  nphi       =  np.random.poisson(phi, len(vs))
  nFv        =  planck * lightspeed * nphi / vs

  xx         =  get_appmags(vs, nFv, {band: filters[band]}, printit=False)

  mags.append(xx[band])

mags = np.array(mags)

pl.hist(mags, bins=50)
pl.axvline(depths[band], c='r')
pl.show()

