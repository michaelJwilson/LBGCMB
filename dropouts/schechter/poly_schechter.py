import numpy as np
import pylab as pl

'''
band   = 'Malkan'
mlim   = 24.6
'''
'''
band   = 'g'
mlim   = 25.8
'''

band   = 'r'
mlim   = 25.8

mm, nn = np.loadtxt('dat/%sDrop.dat' % band, unpack=True)

nn     = nn[mm > 22.5]
mm     = mm[mm > 22.5]

zz     = np.polyfit(mm, np.log10(nn), 4)
p      = np.poly1d(zz)

pl.semilogy(mm, nn)
pl.semilogy(mm, 10. ** p(mm), 'r--')
pl.show()

print(band, 10. ** p(mlim))
