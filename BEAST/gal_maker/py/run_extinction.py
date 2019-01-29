import pylab as pl
import numpy as np

from extinction import calzetti00
from extinction import apply as ext_apply


wave = np.arange(9.e2, 2.5e4)
Fv   = np.ones_like(wave)

##  AV:  extinction in magnitudes at characteristic V band wavelength.
##  RV:  Ratio of total to selective extinction, A_V / E(B-V)

##  pl.plot(wave, Fv, 'k-')

##  'Redden' flux by A_V = 1.0
##  Fv   = ext_apply(calzetti00(wave, a_v=1.0, r_v=4.05, unit='aa'), Fv)
##  pl.plot(wave, Fv, 'c-')

Al  = calzetti00(wave, a_v=1.0, r_v=4.05, unit='aa')  ##  Extinction [mags].
kl  = ext * 4.05 / 1.0                                ##  Alambda = k * Av / Rv;  See http://webast.ast.obs-mip.fr/hyperz/  
                                                      ##  i.e.  E(B-V) = AV / RV; -0.1 to 0.5 
pl.semilogx(wave, kl)

pl.show()
