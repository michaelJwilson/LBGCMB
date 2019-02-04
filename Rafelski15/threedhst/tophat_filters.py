import numpy          as   np
import pylab          as   pl
import astropy.units  as   u

from   collections import  OrderedDict


def tophat_filters(plotit=False):
  filters  = OrderedDict()

  width    = 100.
  
  owave    = np.arange(3.e3, 1.e4,  0.1)
  centers  = np.arange(3.e3, 1.e4, 1.e3)
  
  for i, center in enumerate(centers):
    Ts    = np.ones_like(owave)
    
    Ts[owave < center - width] = 0.0
    Ts[owave > center + width] = 0.0 

    filters[str(i)]       =  {}
    
    temp                  =  owave * u.AA

    filters[str(i)]['ls'] =  owave
    filters[str(i)]['vs'] =  temp.to(u.Hz, equivalencies=u.spectral()).value 
    filters[str(i)]['Ts'] =  Ts

    temp                  =  center * u.AA
    filters[str(i)]['Cs'] =  temp.to(u.Hz, equivalencies=u.spectral()).value

    if plotit:
      pl.plot(owave, Ts, label=str(i))
      pl.show()

  return  filters
