import  matplotlib         as      mpl
import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    pylab              import  rcParams
from    utils              import  latexify


latexify(columns=2, equal=False, fontsize=10, ratio=None, ggplot=True, usetex=True)

LBGCOLORS = ['u-g', 'g-r', 'r-i', 'i-z']

for j, (depth, alpha) in enumerate(zip(['', '_mag_24p5'], [1.0, 0.4])):
 ldepth         = depth
 
 for i, (fname, linestyle) in enumerate(zip(['uvudf', 'aegis'], ['-', '--'])):
  data          = np.loadtxt('../dat/' + fname + '_mag_24p5' + '_colors.txt')
  
  if depth == '_mag_24p5':
   zero_col     = 4

  else:
   zero_col     = 0
   
  idx           =        data[:,0].astype(np.int)
  zee           =        data[:,1]

  umg           =        data[:, 2 + zero_col + 0]
  gmr           =        data[:, 2 + zero_col + 1]
  rmi           =        data[:, 2 + zero_col + 2]
  imz           =        data[:, 2 + zero_col + 3]

  ##  GOLDRUSH z=4 dropouts.
  grush_colorcut  = (gmr > 1.0) & (rmi < 1.0) & (gmr > 1.5 * rmi + 0.8)

  ##  CARS z=3 dropouts.
  cars_colorcut   = (umg > 1.5) & (gmr > -1.0) & (gmr < 1.2) & (umg > 1.5 * gmr + 0.75)

  for colorcut, label, color in zip([grush_colorcut, cars_colorcut], [r'$g$', r'$u$'], ['g', 'b']):
    iidx          = idx[colorcut]  
    zzee          = zee[colorcut]
  
    Ntot          = len(iidx)

    ## Histogram of color cut redshifts. 
    (dNdz, bins)  = np.histogram(zzee, bins = np.arange(0.0, 6.1, 0.10))

    dNdz          = dNdz.astype(np.float)
    bins          = bins[:-1]

    dz            = bins[1] - bins[0]
    midz          = bins + dz / 2.

    join          = ''

    if depth == '_mag_24p5':
      styles      = ['^', 's']
      ldepth      = ' 24.5'

      join        = '_'

    outname       = label + '-' + fname.upper() + join + ldepth
    outname       = '../dat/dNdz/' + outname.translate(None, ''.join(['$', "'", '-', ' '])) + '.txt'

    label         = label + '-' + fname.upper() + ldepth

    pl.plot(midz, dNdz, linestyle, label = label + ': ' + r'%d / %d' % (dNdz.sum(), Ntot), markersize=3, color=color, alpha = alpha)

    np.savetxt(outname, np.c_[midz, dNdz / Ntot / dz], header='  z  dN/dz', fmt='%.6le')

pl.xlabel(r'$z$')
pl.ylabel(r'$dN/dz$')

pl.xlim(  0.0, 6.000)
pl.ylim(-10.0, 6.0e2)

pl.legend(ncol=2, loc=2)

pl.savefig('../plots/dNdz.pdf', bbox_inches='tight')

