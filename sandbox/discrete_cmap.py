import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as pl
import numpy as np
import matplotlib.colors as colors

from mpl_toolkits.axes_grid1 import make_axes_locatable

x       = np.linspace(-2,2,500)
y       = np.linspace(-2,2,500)
XX, YY  = np.meshgrid(x, y)
Z       = np.sin(XX) * np.cos(YY)

cmap        = plt.cm.tab20c ##_r
cmap        = colors.ListedColormap([cmap(i) for i in np.array([2, 1, 4, 5, 8, 9, 12, 13, 16, 17])[::-1]])


bounds      = np.arange(1, 10, 1)
norm        = colors.BoundaryNorm(bounds, cmap.N, clip=False)

plt.pcolormesh(x,y,Z, cmap=cmap, norm=norm)
## plt.colorbar()


ax      = pl.gca()
divider = make_axes_locatable(ax)
cax     = divider.append_axes("right", size="5%", pad=0.05)

cb  = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds+0.5,\
                                boundaries=bounds, format='%.2lf')

cb.set_ticklabels(bounds, update_ticks=True)

plt.show()
