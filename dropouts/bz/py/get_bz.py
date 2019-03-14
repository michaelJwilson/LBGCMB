import  os
import  numpy              as  np
import  pylab              as  pl
import  matplotlib.pyplot  as  plt
import  numpy.ma           as  ma
import  pylab              as  pl
import  pandas             as  pd

from    scipy.interpolate  import  UnivariateSpline, interp1d
from    utils              import  latexify
from    scipy.optimize     import  curve_fit
from    growth_rate        import  growth_factor
from    uncertainties      import  ufloat, unumpy


def bz_fitmodel(z, A, B):
  aa = 1. / (1. + z)
  DD = growth_factor(aa)

  return  A / DD + B / DD ** 2. 

@np.vectorize
def bz_callmodel(z, mlim):
  params = np.array([[24.0, 1.494, 0.077], [24.5, 0.668, 0.195], [25.0, -0.069, 0.320], [25.5, -0.365, 0.356]])
  indexs = np.abs(params[:,0] - mlim) == np.abs(params[:,0] - mlim).min()

  A      = params[indexs, 1]
  B      = params[indexs, 2]

  return  bz_fitmodel(z, A, B)

class drop_bz():
  def __init__(self, type='harikane4'):
    root               = os.environ['LBGCMB'] + '/dropouts/bz/dat/'

    har                = np.loadtxt(root + '%s.dat' % type)
    hare               = np.loadtxt(root + '%s_err.dat' % type)
    
    zs                 = har[:,0][1:]
    mags               = har[0,:][1:]
    
    har[har < -98.9]   = np.NaN
    hare[hare < -98.9] = np.NaN

    har                = ma.masked_invalid(har)[1:,  1:]
    hare               = ma.masked_invalid(hare)[1:, 1:]

    cut                = len(mags)

    harl               = hare[:,:cut]
    haru               = hare[:,cut:]

    self.zs            = zs
    self.ms            = mags 

    self.bs            = har

    self.el            = harl
    self.eu            = haru

    self.zlo           = self.zs.min()
    self.zhi           = self.zs.max()

    self.mlo           = self.ms.min()
    self.mhi           = self.ms.max()

  def print(self):
    print(self.zs)
    print(self.ms)
    print(self.bs)
    print(self.el)
    print(self.eu)
    print(self.zlo)
    print(self.zhi)
    print(self.mlo)
    print(self.mhi)

  def plot(self, marker='o', labelit=True, show=False):
    colors  = plt.rcParams['axes.prop_cycle'].by_key()['color'] 

    for i, m in enumerate(self.ms):
      yerr  = np.array([[self.el[:,i], self.eu[:,i]]])[0,:,:]
      yerr  = np.abs(yerr)

      nudge = np.random.uniform(-1., 1., len(self.zs)) * 1.e-1

      if labelit:
        label = str(m)
      
      else:
        label = ''

      pl.errorbar(self.zs + nudge, self.bs[:,i], yerr=yerr, fmt=marker, label=label, markersize=4, c=colors[i])

    if show:
      pl.legend(frameon=False, ncol=3, loc=2)
      pl.show()

  def get_bzm(self, z, m, printit=False):
    if (z >= self.zlo) & (z <= self.zhi) & (m >= self.mlo) & (m <= self.mhi):
      zindex = np.abs(self.zs - z) == np.abs(self.zs - z).min()
      mindex = np.abs(self.ms - m) == np.abs(self.ms - m).min()

      if printit:
        print(self.bs[zindex, mindex], self.el[zindex, mindex], self.eu[zindex, mindex])

      return  self.bs[zindex, mindex], self.el[zindex, mindex], self.eu[zindex, mindex]
    
    else:
      raise UserWarning('\n\nHarikane bias not available at z=%.2lf' % z)

  def join(self, x):
    result     = drop_bz()

    result.zs  = np.concatenate([self.zs, x.zs])

    assert  len(result.zs) == len(np.unique(result.zs))
    assert  np.allclose(self.ms, x.ms)

    ##  Sort.
    indices    = np.argsort(result.zs, kind='quicksort')

    result.zs  = result.zs[indices]    
    
    result.zlo = result.zs.min()
    result.zhi = result.zs.max()

    result.mlo = result.ms.min()
    result.mhi = result.ms.max()

    result.bs  = np.vstack((self.bs, x.bs))
    result.bs  = result.bs[indices, :]

    result.el  = np.vstack((self.el, x.el))
    result.el  = result.el[indices, :]

    result.eu  = np.vstack((self.eu, x.eu))
    result.eu  = result.eu[indices, :]

    return  result

  def mask_mlim(self, mlim):
    indices    = (np.abs(self.ms - mlim) != np.abs(self.ms - mlim).min())

    self.ms    = self.ms[indices]

    self.bs    = self.bs[:, indices]

    self.el    = self.el[:, indices]
    self.eu    = self.eu[:, indices]

    self.mlo   = self.ms.min()
    self.mhi   = self.ms.max()

  def mask_z(self, zee):
    indices    = (np.abs(self.zs - zee) != np.abs(self.zs - zee).min())

    self.zs    = self.zs[indices]

    self.bs    = self.bs[indices, :]

    self.el    = self.el[indices, :]
    self.eu    = self.eu[indices, :]

    self.zlo   = self.zs.min()
    self.zhi   = self.zs.max()

  def fit(self, mlim, plotit=False):
    mindex     = np.abs(self.ms - mlim) == np.abs(self.ms - mlim).min()

    self.bs    = np.ma.masked_invalid(self.bs)
    
    xx         = self.zs[np.array(self.bs[:, mindex].nonzero()).T[:,0]]
    yy         = self.bs[:, mindex][self.bs[:, mindex].nonzero()]
    ee         = np.abs(self.el[:, mindex][self.bs[:, mindex].nonzero()])

    popt, pcov = curve_fit(bz_fitmodel, xx, yy, p0 = None, sigma = ee, absolute_sigma = False)

    if plotit:
      colors   = plt.rcParams['axes.prop_cycle'].by_key()['color']

      xxs      = np.arange(2.5, 6.0, 0.01)
      pl.plot(xxs, bz_fitmodel(xxs, *popt), c=colors[np.where(mindex == True)[0][0]])

    return  popt, pcov


if __name__ == '__main__':
    print('\n\nWelcome.\n\n')

    x = drop_bz(type='cars3')  ##  ['cars3', 'harikane4']
    y = drop_bz(type='harikane4')

    z = x.join(y)

    z.mask_mlim(26.0)
    z.mask_mlim(26.5)

    z.mask_z(5.9)    

    for mlim in [24.0, 24.5, 25.0, 25.5]:
      popt, pcov = z.fit(mlim, plotit=True)

    pl.xlim(2.5, 5.5)
    pl.ylim(2.0, 9.0)

    z.plot(show=False, marker='s')

    ##  x.plot(show=False, marker='s')
    ##  y.plot(show=True,  marker='o')

    zs = np.arange(2.5, 5.5, 0.01)
    bs = bz_callmodel(z=zs, mlim=24.)
    
    pl.plot(zs, bs, 'c-')

    pl.legend()
    pl.show()
    
    print('\n\nDone.\n\n')
