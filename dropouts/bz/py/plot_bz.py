import  os
import  numpy              as  np
import  pylab              as  pl
import  matplotlib.pyplot  as  plt

from    scipy.interpolate  import  UnivariateSpline, interp1d
from    utils              import  latexify
from    scipy.optimize     import  curve_fit
from    growth_rate        import  growth_factor
from    get_bz             import  drop_bz

latexify(fig_height=2., fig_width=3., columns=2, equal=False, fontsize=10, ratio=0.5, ggplot=True, usetex=True)

##  ii   = interp1d(har[:-1,0], har[:-1,1], 'linear', copy=True, bounds_error=False, fill_value='extrapolate', assume_sorted=False)
##  ii  = UnivariateSpline(har[:-1,0], har[:-1,1], k=3, s=2, ext=0, check_finite=False)

x = drop_bz(type='cars3')
y = drop_bz(type='harikane4')
z = x.join(y)

z.mask_mlim(26.0)
z.mask_mlim(26.5)
z.mask_z(5.9)

x.mask_mlim(26.0)
x.mask_mlim(26.5)
x.mask_z(5.9)

y.mask_mlim(26.0)
y.mask_mlim(26.5)

result = []

for mlim in [24.0, 24.5, 25.0, 25.5]:
  popt, pcov = z.fit(mlim, plotit=True)

  result.append([mlim, *popt])

result = np.array(result)

np.savetxt('../dat/bz_modelparams.dat', result, fmt='%.3lf')

x.plot('o', labelit=False)
y.plot('s', labelit=True, show=False)

##  z.plot(show=False, marker='s')

pl.xlim(2.75, 5.10)
pl.ylim(1.75, 9.75)

pl.xlabel(r'$z$')
pl.ylabel(r'$b(z)$')

pl.legend(loc=2, ncol=2, handletextpad=.05, frameon=False)

plt.tight_layout()

##  pl.show()
pl.savefig(os.environ['LBGCMB'] + '/dropouts/bz/plots/bz.pdf') 
