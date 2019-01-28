import  pickle
import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    numpy.linalg       import  inv
from    scipy.interpolate  import  interp1d
from    params             import  get_params
from    cosmo              import  cosmo
from    utils              import  prefactor, comoving_distance
from    Cgg                import  Ngg
from    scipy.misc         import  derivative
from    growth_rate        import  growth_factor
from    pmh                import  linz_bz
from    qso_bz             import  qso_bz
from    Gaussian_pz        import  Gaussian_pz
from    scipy.special      import  erf
from    prep_pk            import  prep_pk
from    pnw                import  pnw
from    fisher_contour     import  plot_ellipse
from    utils              import  latexify, sci_notation


def Psp(z, b1, b2):
  ks, Ps = prep_pk(z)
  
  return  ks, b1(z) * b2(z) * Ps

def pz_zp(z, sigma, zlo, zhi):
  ##  prob(z | zp) for zp \in [zlo, zhi] defined by spec. slice
  ##  See eqn. (11) of 1302.6015
  sigma *= (1. + np.mean([zlo, zhi]))

  xhi    = (zhi - z) / np.sqrt(2.) / sigma
  xlo    = (zlo - z) / np.sqrt(2.) / sigma

  return  0.5 * (erf(xhi) - erf(xlo))

@np.vectorize
def pzp(z, z0=0.4):
  ##  prob(zp)
  ##  See eqn. (12) of 1302.6015   
  return  z * z / 2. / z0 / z0 / z0 * np.exp(-z / z0)

@np.vectorize
def pz_photo(z, sigma, z0, zlo, zhi):
  dzp  = 0.01
  zp   = np.arange(0.0, 5.0, dzp)

  return  pzp(z, z0=0.4) * pz_zp(z, sigma, zlo, zhi)

def Cpp():
  return  None

def Css(zlo, zhi, b1, noks=False):
  ##  In the zeff approx.
  lochi    = comoving_distance(zlo)
  hichi    = comoving_distance(zhi)

  midz     = np.mean([zlo, zhi])
  ks, Ps   = Psp(midz, b1, b1)

  ##  Fraction of spectroscopic galaxies that reside in the 
  ##  bin is unity by definition. 
  Cs       = Ps / (hichi - lochi) 

  if noks:
    return Ps
  
  else:
    return ks, Ps 

@np.vectorize
def Csp(z0, sigma, zlo, zhi, b1, b2, alpha=1.0, disp = 4.0, nowiggle=False, nobroadband=True, noks=False):
  ##  In the zeff approx.
  lochi    = comoving_distance(zlo)
  hichi    = comoving_distance(zhi)

  midz     = np.mean([zlo, zhi])
  
  if nowiggle == 'damped':
    ##  BAO work.
    ks, Ps = Psp(midz, b1, b2)
    ks, Ns = pnw(midz, b1, b2)

    Ps     = interp1d(ks, Ps, kind='linear', copy=True, bounds_error=False, fill_value=0.0)
    Ns     = interp1d(ks, Ns, kind='linear', copy=True, bounds_error=False, fill_value=0.0)

    kp     = np.copy(ks)  ##  Wavenumber in fiducial cosmology.
    ks     = alpha * kp   ##  Wavenumber in true cosmology, kp = (k / alpha).

    if nobroadband:
      Ps     = (Ps(kp) - Ns(kp)) * np.exp(- ks * ks * disp * disp / 2.) + Ns(ks)
      
    else:
      ##  Argument change to broadband term. 
      Ps     = (Ps(kp) - Ns(kp)) * np.exp(- ks * ks * disp * disp / 2.) + Ns(kp)

    ##  Rename to match rest of script.
    ks     = kp

  elif nowiggle:
    ks, Ps = pnw(midz, b1, b2) 

  else:
    ##  Full wiggle.
    ks, Ps = Psp(midz, b1, b2)

  dz      = 0.01
  zs      = np.arange(0.0, 10.0, dz)

  nbar    = dz * np.sum(pz_photo(zs, sigma, z0, zlo, zhi))
  
  zs      = np.arange(zlo, zhi + dz, dz)
  frac    = dz * np.sum(pz_photo(zs, sigma, z0, zlo, zhi))

  Cs      = frac * Ps / (hichi - lochi)

  if noks:
    return  Cs

  else:
    return  ks, Cs

def deriv_Csp(param, z0, sigma, zlo, zhi, b1, b2, alpha=1.0, disp=4.0, nowiggle=False, nobroadband=True):
  ##  Derivative of Csp wrt alpha or displacement (Sigma).
  if param == 'alpha':
    fz   = lambda z:  Csp(z0, sigma, zlo, zhi, b1, b2, alpha=z, disp = disp, nowiggle='damped', nobroadband=nobroadband, noks=True)
    dfdz = derivative(fz, alpha, dx=0.001, n=1, order=3)

  elif param == 'disp':
    fz   = lambda z:  Csp(z0, sigma, zlo, zhi, b1, b2, alpha=alpha, disp = z, nowiggle='damped', nobroadband=nobroadband, noks=True)
    dfdz = derivative(fz, disp, dx=0.1, n=1, order=3)

  else:
    raise UserWarning('\n\nRequested derivative is not available.\n\n')

  return dfdz

def cov_Csp(z0, sigma, zlo, zhi, b1, b2, ns):
  ##  Currently assume dominated by cross-correlation term.
  ##  See eqn. (14) of 1302.6015
  ks, Cxy =  Csp(z0, sigma, zlo, zhi, b1, b2, nowiggle=True)

  ks, Cxx =  Css(zlo, zhi, b1, noks=False)
  Cyy     =  np.zeros_like(Cxx)

  lochi   =  comoving_distance(zlo)
  hichi   =  comoving_distance(zhi)
  
  pdist   =  hichi - lochi
  ns     *=  pdist

  ##  \bar n_p e zs;  See eqn.(14).
  dz      = 0.01
  zs      = np.arange(0.0, 10.0, dz)

  _np     = dz * np.sum(pz_photo(zs, sigma, z0, zlo, zhi))

  return  Cxy ** 2. + (Cxx + 1. / ns) * (Cyy + 1. / _np) 

def nmodes2D(ks, dk, zlo, fsky):
  return  2. * fsky * ks * dk * (1. + zlo) ** 2. * comoving_distance(zlo) ** 2. 

def TDdChi2(V, zspec, kmax, nbar, b1, dk=0.01):
  ks, Ps    =  Psp(zspec, b1, b2)
  ks, PNs   =  pnw(zspec, b1, b2)

  Ps        =  Ps[ks <= kmax]
  PNs       = PNs[ks <= kmax]
  ks        =  ks[ks <= kmax]

  result    = (1. / (2. * np.pi) ** 2.) * ks * ks * dk * ((Ps - PNs) / (Ps + 1. / nbar)) ** 2.
  result    = np.sum(result)
  result   *= V               ##  eqn. (18) of 1302.6015 has wrong V factor for mode counting?

  return  result


if __name__ == '__main__':
  print("\n\nWelcome to Nishizwa.\n\n")

  fsky    =  0.1
  
  zspec   =  0.8
 
  zlo     =  0.6
  zhi     =  1.0

  b1      = qso_bz
  b2      = linz_bz

  z0      =  0.4
  sigma   =  1.0
  
  '''
  ks, Ps  = Psp(zspec, b1, b2)
  '''
  '''
  result  = []
  kmaxs   = np.arange(0.1, 0.5, 0.05)

  for kmax in kmaxs:
    dChi2 = TDdChi2(1.e9, zspec, kmax=kmax, nbar=1.e-3, b1=b1)

    print('%.2lf \t %.2lf' % (kmax, dChi2))

    result.append(dChi2)

  result  = np.array(result)  

  pl.plot(kmaxs, result)
  pl.show()
  '''
  '''
  p1      = lambda z:  pp(z, sigma, 1.0, 1.8) 
  p2      = lambda z:  Gaussian_pz(z, z0 = 2.96, sigma = 0.24)

  zs      = np.arange(0.0, 5.0, 0.01)

  pl.plot(zs, pzp(zs, z0=z0), 'k-', label=r'$Prob(z_p)$')

  for [x,y] in [[0.6, 1.0],[1.0, 1.8],[1.8, 3.2]]:
    pl.plot(zs, pz_photo(zs, sigma, z0, x, y), label='1')
    pl.plot(zs, pz_photo(zs, sigma, z0, x, y), label='2')
    pl.plot(zs, pz_photo(zs, sigma, z0, x, y), label='3')

  pl.legend()
  pl.show()
  '''
  '''
  ks, Cs  = Csp(z0, sigma, zlo, zhi, b1, b2)
  ls      = ks * comoving_distance(1.4)

  pl.loglog(ls, prefactor(ls, n=2) * Cs)
  pl.show()
  '''
  '''
  result = []
  sigmas = [100., 0.1]

  for sigma in sigmas:
    ks, Cs  = Csp(z0, sigma, zlo, zhi, b1, b2, nowiggle=False)
    ks, Ns  = Csp(z0, sigma, zlo, zhi, b1, b2, nowiggle=True)

    dChi2   = np.sum(nmodes2D(ks, 0.01, zlo, fsky) * (Cs - Ns) ** 2. / cov_Csp(z0, sigma, zlo, zhi, b1, b2))

    pl.plot(ks, nmodes2D(ks, 0.01, zlo, fsky) * (Cs - Ns) ** 2. / cov_Csp(z0, sigma, zlo, zhi, b1, b2), label=r'$%.1lf, \ %.1lf$' % (sigma, dChi2))

    result.append(dChi2)

  pl.legend()
  pl.show()

  result = np.array(result)
  '''
  ##  pl.plot(sigmas, result)
  ##  pl.show()

  Fish   = np.zeros(4).reshape(2,2)

  for i, x in enumerate(['alpha', 'disp']):
    for j, y in enumerate(['alpha', 'disp']):
      d1          = deriv_Csp(x, z0, sigma, zlo, zhi, b1, b2, alpha=1.0, disp=4.0, nowiggle='damped', nobroadband=True)
      d2          = deriv_Csp(y, z0, sigma, zlo, zhi, b1, b2, alpha=1.0, disp=4.0, nowiggle='damped', nobroadband=True)

      ##  Contribution of one redshift slice; ns = comoving number density of spec. z galaxies. 
      Fish[i,j]   = np.sum(d1 * d2 / cov_Csp(z0, sigma, zlo, zhi, b1, b2, ns=1.e-4))  

  iFish     = inv(Fish)

  sig_alpha = np.sqrt(iFish[0,0])
  sig_disp  = np.sqrt(iFish[1,1])

  print('\n\nSolved for: sig. alpha:  %.4le \t sig. disp.:  %.4le.\n\n' % (sig_alpha, sig_disp))

  ##  And plot contour ...                                                                                                                                 
  pl.clf()

  ##  latexify(columns=1, equal=True)

  fig    = plt.gcf()
  ax     = plt.gca()

  for mass_level, color, alpha in zip([0.99, 0.95, 0.68], ['b', 'b', 'b'], [0.2, 0.4, 0.6]):
    plot_ellipse(x_cent = 1.0, y_cent = 4.0, ax = ax, cov = iFish, mass_level = mass_level,\
                 fill=True, fill_kwargs={'alpha': alpha, 'c': color}, plot_kwargs={'c': color, 'alpha': 0.0})

  str    = sci_notation(sig_alpha, decimal_digits=2, precision=None, exponent=None)

  pl.plot(1.0, 4.0, 'k*', markersize=5, label=r'$(\sigma_\alpha, \sigma_\Sigma)=$' + '(%.3lf, ' % sig_alpha + '%.0lf)' % sig_disp)

  pl.xlabel(r'$\alpha$')
  pl.ylabel(r'$\Sigma \ [(h^{-1} \rm{Mpc})]$')

  pl.ylim(-400., 400.)

  pl.legend(loc=3)

  plt.tight_layout()

  pl.savefig('plots/nishikawa.pdf')
  
  print("\n\nDone\n\n")
