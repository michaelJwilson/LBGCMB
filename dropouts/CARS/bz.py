import numpy as np
import pylab as pl


##  GOLDRUSH
dat = np.loadtxt('harikane.dat')

pl.scatter(dat[:,0], dat[:,1], c='b',      marker='^', label='24.0', s=7)
pl.scatter(dat[:,0], dat[:,2], c='orange', marker='^', label='', s=7)
pl.scatter(dat[:,0], dat[:,3], c='c',      marker='^', label='', s=7)
pl.scatter(dat[:,0], dat[:,4], c='y',      marker='^', label='', s=7)

##  CARS                                                                                                
dat = np.loadtxt('cars.dat')

pl.scatter(dat[:,0], dat[:,1], c='orange', marker='s', label='24.5', s=7)
pl.scatter(dat[:,0], dat[:,2], c='c',      marker='s', label='25.0', s=7)
pl.scatter(dat[:,0], dat[:,3], c='y',      marker='s', label='25.5', s=7)
pl.scatter(dat[:,0], dat[:,4], c='k',      marker='s', label='26.0', s=7)
pl.scatter(dat[:,0], dat[:,5], c='g',      marker='s', label='26.5', s=7)


pl.ylim(1.75, 9.0)

pl.xlabel(r'$z$')
pl.ylabel(r'$b(z)$')

pl.legend(loc=2, ncol=3)

pl.savefig('bz.pdf') 
