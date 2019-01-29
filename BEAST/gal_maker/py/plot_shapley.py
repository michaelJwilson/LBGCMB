import os
import numpy as np
import pylab as pl

from   lbg_maker    import lbg_maker
from   prep_filters import prep_filters
from   madau        import lephare_madau
from   extinction   import calzetti00
from   extinction   import apply as ext_apply


root              =  os.environ['BEAST']

redshifts         =   1.5 + np.linspace(0., 4.0, 4)
magnitudes        =  23.5 * np.ones_like(redshifts)

flux, wave, meta  =  lbg_maker(ngal=None, restframe=False, printit=False, test=True, seed=314, redshifts=None, magnitudes=None, hiwave=2.e4)

filters           =  prep_filters(['LSST', 'VIDEO'], normed=True)

for i, band in enumerate(list(filters.keys())):
    pl.fill_between(filters[band]['ls'], filters[band]['Ts'] * 1., alpha=0.3, label='')  ##  label=filters[band]['ppkey']

for i, x in enumerate(flux):
    rwave = wave / (1. + meta['REDSHIFT'][i])
    x     = ext_apply(calzetti00(rwave, a_v=0.2, r_v=4.05, unit='aa'), x)

    ## Madau extinction. 
    x    *= lephare_madau(rwave, meta['REDSHIFT'][i]) 

    pl.semilogy(wave, x, label=r'$%+ .1lf\AA, z=%.1lf$' % (meta['Lya-EW'][i], meta['REDSHIFT'][i]))

pl.xlabel(r'$\lambda \ [\rm{Angstroms}]$')
pl.ylabel(r'$F_\lambda \ [10^{-17} \ \rm{ergs} / s / \rm{cm}^2 / \AA]$')

pl.xlim(3.e3, 1.9e4)
pl.ylim(0.1,  100.1)

pl.legend(ncol=1, loc=1)
pl.savefig(root + '/gal_maker/plots/shapley.pdf')
