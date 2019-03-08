import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

from   scipy.interpolate import UnivariateSpline, interp1d
from   utils             import latexify


latexify(fig_height=2., fig_width=3., columns=2, equal=False, fontsize=10, ratio=0.5, ggplot=True, usetex=True)

##  GOLDRUSH
har = np.loadtxt('harikane.dat')

zs  = np.arange(2.7, 5.1, 0.01) 
ii  = interp1d(har[:-1,0], har[:-1,1], 'linear', copy=True, bounds_error=False, fill_value='extrapolate', assume_sorted=False)
##  ii  = UnivariateSpline(har[:-1,0], har[:-1,1], k=3, s=2, ext=0, check_finite=False)

pl.plot(har[:-1,0], har[:-1,1], '-', c='b', label='')
pl.plot(zs, ii(zs), '--', c='b', label='')

##  CARS                                                                                                
cars     = np.loadtxt('cars.dat')

##  Orange
ozs, obs = np.concatenate([cars[:1,0], har[:-1,0]]), np.concatenate([cars[:1,1], har[:-1,2]]) 
ii       = interp1d(ozs, obs, 'quadratic', copy=True, bounds_error=False, fill_value='extrapolate', assume_sorted=False)
##  ii   = UnivariateSpline(ozs, obs, k=3, s=2, ext=0, check_finite=False)

pl.plot(zs, ii(zs), 'orange')

##  Cyan.
ozs, obs = np.concatenate([[cars[0,0]], [cars[1,0]], [har[0,0]], har[1:,0]]), np.concatenate([[cars[0,2]], [cars[1,2]], [har[0,3]], har[1:,3]])
##  ii   = interp1d(ozs, obs, 'quadratic', copy=True, bounds_error=False, fill_value='extrapolate', assume_sorted=False)
ii       = UnivariateSpline(ozs, obs, k=3, s=2, ext=0, check_finite=False)

pl.plot(zs, ii(zs), 'cyan')

##  Yellow.                                                                                                                                    
ozs, obs = np.concatenate([[cars[0,0]], [cars[1,0]], [har[0,0]], har[1:-1,0]]),\
           np.concatenate([[cars[0,3]], [cars[1,3]], [har[0,4]], har[1:-1,4]])

##  ii   = interp1d(ozs, obs, 'quadratic', copy=True, bounds_error=False, fill_value='extrapolate', assume_sorted=False)
ii       = UnivariateSpline(ozs, obs, k=2, s=5, ext=0, check_finite=False)

pl.plot(zs, ii(zs), 'y')

##
pl.scatter(har[:,0], har[:,1], c='b',      marker='^', label='24.0', s=14)
pl.scatter(har[:,0], har[:,2], c='orange', marker='^', label='',     s=14)
pl.scatter(har[:,0], har[:,3], c='c',      marker='^', label='',     s=14)
pl.scatter(har[:,0], har[:,4], c='y',      marker='^', label='',     s=14)

##
pl.scatter(cars[:,0], cars[:,1],     c='orange', marker='s', label='24.5', s=14)
pl.scatter(cars[:,0], cars[:,2],     c='c',      marker='s', label='25.0', s=14)
pl.scatter(cars[:-1,0], cars[:-1,3], c='y',      marker='s', label='25.5', s=14)

pl.xlim(2.90, 5.10)
pl.ylim(1.75, 9.75)

pl.xlabel(r'$z$')
pl.ylabel(r'$b(z)$')

pl.legend(loc=2, ncol=2, handletextpad=.05)

plt.tight_layout()

##  pl.show()
pl.savefig('bz.pdf') 
