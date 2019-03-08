import numpy             as np
import pylab             as pl
import matplotlib.pyplot as plt


x = np.arange(10)

pl.plot(x, x, c='w')

frame1 = plt.gca()

frame1.axes.xaxis.set_ticklabels([])
frame1.axes.yaxis.set_ticklabels([])

plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,        # ticks along the bottom edge are off
    right=False,       # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off

pl.savefig('plots/blank.pdf')
