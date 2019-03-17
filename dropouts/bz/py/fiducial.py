import numpy as np
import pylab as pl

from   get_bz import bz_callmodel


##  Photometric.
##  print(bz_callmodel(2.2, 24.0))
print(bz_callmodel(3.0, 24.6))
print(bz_callmodel(4.0, 25.8))
print(bz_callmodel(5.0, 25.8))

print('\n\n')

##  Spectroscopic.
print(bz_callmodel(2.2, 24.0))
print(bz_callmodel(3.0, 24.0))
print(bz_callmodel(4.0, 25.5))
print(bz_callmodel(5.0, 25.5))
