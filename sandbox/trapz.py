import numpy as np
import pylab as pl


xlo =    1.0
xhi =  100.0

xs  =    np.arange(xlo, xhi, 0.01)
xs  =    np.arange(np.sqrt(xlo), np.sqrt(xhi), 0.01) ** 2.
xs  =  np.logspace(np.log10(xlo), np.log10(xhi), 1000)

ys  = xs ** 2.

I   = np.trapz(ys, xs)
J   = (xhi ** 3. - xlo ** 3.) / 3.

print(I, J)
