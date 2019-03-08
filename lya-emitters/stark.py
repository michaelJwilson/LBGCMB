import numpy             as np 
import pylab             as pl
import astropy.units     as u
import matplotlib.pyplot as plt

from   astropy.cosmology import FlatLambdaCDM
from   utils             import latexify


h     = 0.67


cosmo = FlatLambdaCDM(H0=100. * h, Om0=0.26, Tcmb0=2.725)


zs    = np.array([3, 4., 5.])
ms    = np.arange(20., 30., 0.01)

latexify(columns=1, equal=True, fontsize=12, ratio=None, ggplot=True, usetex=True)

for z in zs:
    DL  = cosmo.luminosity_distance(z).value * 1.e6  ## parsec.
    Ms  = ms + 2.5 * np.log10(1. + z) - 5. * np.log10(DL / 10.)

    x   = 1.09 + 0.047 * Ms

    lim = -20.51

    x[Ms > lim]  = 5.46 + 0.260 * Ms[Ms > lim]

    pl.plot(ms, x, label=r'$z=%.0lf$' % z)

pl.xlim(24.,  27.)
pl.ylim(0.,   0.8)

pl.legend(loc=2, frameon=False)

pl.xlabel(r'$m_{UV}$')
pl.ylabel(r'$X_{\alpha}$')

plt.tight_layout()

##  pl.show()
pl.savefig('plots/frac.pdf')
