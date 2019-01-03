import numpy as np
import pylab as pl


mAB = np.arange(-24., 24., 0.5)

Fv  = - (mAB + 48.6) / 2.5 
Fv  = 10. ** Fv

pl.semilogy(mAB, Fv)

pl.xlabel(r'$m_{AB}$')
pl.ylabel(r'$\log_{10} \ | \ F_{\nu} \ / \ \rm{erg} \ / \ s \ / \ \rm{cm}^2 \ / \ \rm{Hz} \ | $')

pl.show()
