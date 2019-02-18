import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

from   scipy.interpolate import UnivariateSpline, interp1d
from   utils             import latexify


latexify(fig_width=None, fig_height=None, columns=1, equal=True, fontsize=10, ratio=None, ggplot=True, usetex=True)

##  GOLDRUSH
har = np.loadtxt('harikane.dat')

zs  = np.arange(3., 5.1, 0.01) 
ii  = interp1d(har[:-1,0], har[:-1,1], 'linear', copy=True, bounds_error=False, fill_value='extrapolate', assume_sorted=False)

pl.plot(har[:-1,0], har[:-1,1], '-', c='b', label='')
pl.plot(zs, ii(zs), '--', c='b', label='')

pl.scatter(har[:,0], har[:,1], c='b',      marker='^', label='24.0', s=9)
pl.scatter(har[:,0], har[:,2], c='orange', marker='^', label='', s=9)
pl.scatter(har[:,0], har[:,3], c='c',      marker='^', label='', s=9)
pl.scatter(har[:,0], har[:,4], c='y',      marker='^', label='', s=9)

##  CARS                                                                                                
cars = np.loadtxt('cars.dat')

pl.scatter(cars[:,0], cars[:,1], c='orange', marker='s', label='24.5', s=9)
pl.scatter(cars[:,0], cars[:,2], c='c',      marker='s', label='25.0', s=9)
pl.scatter(cars[:-1,0], cars[:-1,3], c='y',      marker='s', label='25.5', s=9)

##  Orange
ozs, obs = np.concatenate([cars[:1,0], har[:-1,0]]), np.concatenate([cars[:1,1], har[:-1,2]]) 
ii       = interp1d(ozs, obs, 'quadratic', copy=True, bounds_error=False, fill_value='extrapolate', assume_sorted=False)

pl.plot(zs, ii(zs), 'orange')

##  Cyan.
ozs, obs = np.concatenate([[cars[0,0]], [har[0,0]], [cars[1,0]], har[1:,0]]), np.concatenate([[cars[0,2]], [har[0,3]], [cars[1,2]], har[1:,3]])
ii       = interp1d(ozs, obs, 'quadratic', copy=True, bounds_error=False, fill_value='extrapolate', assume_sorted=False)
ii       = UnivariateSpline(ozs, obs, k=3, s=2, ext=0, check_finite=False)

pl.plot(zs, ii(zs), 'cyan')

##  Yellow.                                                                                                                                    
ozs, obs = np.concatenate([[cars[0,0]], [har[0,0]], [cars[1,0]], har[1:-1,0]]),\
           np.concatenate([[cars[0,3]], [har[0,4]], [cars[1,3]], har[1:-1,4]])

ii       = interp1d(ozs, obs, 'quadratic', copy=True, bounds_error=False, fill_value='extrapolate', assume_sorted=False)
ii       = UnivariateSpline(ozs, obs, k=2, s=5, ext=0, check_finite=False)

print(ozs)
print(obs)

pl.plot(zs, ii(zs), 'yellow')

pl.xlim(2.90, 5.1)
pl.ylim(1.75, 9.0)

pl.xlabel(r'$z$')
pl.ylabel(r'$b(z)$')

pl.legend(loc=4, ncol=1, handletextpad=.05)
plt.tight_layout()
pl.savefig('bz.pdf') 
