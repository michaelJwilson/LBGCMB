import  os
import  pickle
import  numpy                            as      np
import  pylab                            as      pl
import  matplotlib.pyplot                as      plt

from    utils                            import  latexify
from    sklearn.gaussian_process         import  GaussianProcessRegressor
from    sklearn.gaussian_process.kernels import  RBF, ConstantKernel as C
from    scipy                            import  sparse


latexify(columns=2, equal=False, fontsize=12, ggplot=True, usetex=True, ratio=0.4)

def get_pz(band, depth, zlimit):
  root          = os.environ['LBGCMB']
  midz, pz, spz = np.loadtxt(root + '/dropouts/interlopers/pz/pz_%s_%s.txt' % (band, depth), unpack=True)

  pz            =    pz[midz <= zlimit]
  spz           =   spz[midz <= zlimit]
  midz          =  midz[midz <= zlimit]

  ##  GP                                                                                                                                                                           
  kernel        = C(1.0, (1e-3, 1e3)) * RBF(.1, (1e-2, 1e0))
  gp            = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=15, alpha=spz**2.)

  ##  Fit to data using Maximum Likelihood Estimation of the parameters                                                                                                            
  gp.fit(midz.reshape(-1, 1), pz)

  ## 
  dz            = 0.01
  zs            = np.arange(0.0, 10., dz).reshape(-1, 1)
  norm          = np.sum(gp.predict(zs, return_std=False)) * dz

  print('\n\nCalculated raw normalisation as %.4le' % norm)

  ngp           = lambda z:  gp.predict(np.array(z).reshape(-1,1), return_std=False)[0] / norm 

  return  midz, pz, spz, gp, ngp


if __name__ == '__main__':
  dz     = 0.25
  colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

  for band in ['u', 'g']:
    pl.clf()

    for ii, depth in enumerate(['Full', 'Degraded']):
      midz, pz, spz, gp, ngp = get_pz(band, depth, 5.5)
    
      pl.errorbar(midz, pz, spz, label=depth, color=colors[ii], fmt='', ls='')
    
      zs        = np.arange(0.0, midz.max(), 0.01).reshape(-1, 1)
      ys        = ngp(zs)
      ##  ys, sigma = gp.predict(zs, return_std=True)

      pl.plot(zs, ys, 'k-', alpha=0.3, color=colors[ii])
      
      ##  plt.fill(np.concatenate([zs, zs[::-1]]), np.concatenate([ys - 1.9600 * sigma, (ys + 1.9600 * sigma)[::-1]]), alpha=.2, fc=colors[ii], ec='None', label='')
      ##  sL    = sparse.triu(gp.L_, k=-1)

    ##  pl.xlim(0.0, 6.00)
    pl.ylim(0.0, 1.25)

    pl.legend(ncol=2, loc=2, frameon=False, title=r'$%s$-dropouts' % band)

    pl.xlabel(r'$z$', fontsize=14)
    pl.ylabel(r'$p(z)$')

    plt.tight_layout()

    ##   
    ax  = pl.gca()
    fig = pl.gcf()
    
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_color('black')
    
    plt.tight_layout()
    pl.show()
    ##  pl.savefig('plots/pz_%s.pdf' % band)

print('\n\nDone.\n\n')
