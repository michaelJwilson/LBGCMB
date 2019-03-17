import numpy as np
import pylab as pl


params = np.array([[24.0,  1.178, 0.044],\
                   [24.5,  0.538, 0.115],\
                   [25.0, -0.032, 0.190],\
                   [25.5, -0.262, 0.212]])
deg    =   1
mshift = -25.

x      = params[:,0] + mshift
y      = params[:,1]
z      = params[:,2] 

pc     = np.polyfit(x, y, deg)
pm     = np.poly1d(pc)

qc     = np.polyfit(x, z, deg)
qm     = np.poly1d(qc)

ms     = np.arange(22.0, 27.5, 0.01) + mshift

pl.plot(x, y,   '^',  markersize=4)
pl.plot(x, z,   '^',  markersize=4)

print(pc)
print(qc)

## 
pl.plot(ms, pm(ms), 'c-')
pl.plot(ms, qm(ms), 'm-')
pl.show()
