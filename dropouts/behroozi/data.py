import numpy as np
import pylab as pl

from   utils import latexify


latexify(fig_width=None, fig_height=None, columns=2, equal=False, fontsize=10, ratio=None, ggplot=True, usetex=True)

##  Observational Stellar mass.
for zee in [4, 5]:
  dat        = np.loadtxt('dat/data-compilation/smf_ms/kslee_z%d.smf' % zee)
  Mstar      = 10. ** dat[:,0]  ## [Msun].
  nbarperdex = 10. ** dat[:,1]  ## [Mpc^-3 dex^{-1}]

  pl.loglog(Mstar, nbarperdex, '^', markersize=3, label=r'$z=%.1lf$' % zee)

pl.legend()
pl.xlabel(r'$M_* \ [M_{\odot}]$')
pl.ylabel(r'$\bar{n} \ [\rm{Mpc}^{-3} \ \rm{dex}^{-1}] $')
pl.savefig('plots/obs_sm.pdf', bbox_inches='tight')

