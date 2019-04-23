import  collections
import  numpy                as      np
import  pylab                as      pl
import  matplotlib.pyplot    as      plt
import  matplotlib           as      mpl

from    cosmo                import  cosmo
from    nbar                 import  dVols
from    qso.get_pz           import  get_pz
from    elg_nz               import  get_pz  as elg_pz
from    completeness         import  get_completeness
from    Gaussian_pz          import  Gaussian
from    lensing              import  lensing_kernel
from    cosmo                import  cosmo
from    params               import  get_params
from    scipy                import  misc, ndimage
from    whitebook_pz         import  whitebook_pz, lsst_chang_pz
from    dndz                 import  getpz_H09
from    astropy.convolution  import  convolve, Box1DKernel, Gaussian1DKernel
from    scipy.interpolate    import  interp1d
from    utils                import  latexify


latexify(fig_width=None, fig_height=None, columns=2, equal=False, fontsize=10, ratio=0.5, usetex=True, ggplot=True)

params = get_params()

if __name__ == "__main__":
  from pylab import rcParams

  print("\n\nWelcome to a LSST p(z) calculator.\n\n")

  ax      = pl.gca()

  dz      = 0.1
  zs      = np.arange(0.0, 1100.0, dz)

  ##  lensing kernel                                                                                                                            
  kernel  = lensing_kernel(zs) / (100. * params['h_100'] * cosmo.efunc(zs))
  norm    = np.sum(kernel * dz) 

  ax.fill_between(zs, np.zeros_like(zs), kernel / kernel.max(), alpha=0.15, color='navy', label=r'$W^{\kappa}$')
   
  ##  and galaxies.. 
  dz      = 0.01
  zmax    =  7.0

  zs      = np.arange(0.0, zmax,  dz)
  ilims   = np.arange(20., 35.,  0.5)

  for k, ilim in enumerate(ilims):
    pl.plot(zs, whitebook_pz(zs, ilim), 'k', alpha= 1.0 * k / len(ilims), lw=0.75)
  
  ## Surveys specified by limiting i magnitude; assumes p(z) follows LSST science book form.                                                               
  colors = ['maroon', 'indianred', 'mediumblue', 'teal']

  ilims  = collections.OrderedDict()
  
  ilims['DES']      = 23.5
  ilims['DES-Deep'] = 24.5
  ilims['LSST-Y1']  = 24.1
  ilims['LSST-Y10'] = 25.3
  
  for k, ilim in enumerate(ilims):
    try:
      pl.plot(zs, lsst_chang_pz(zs, ilims[ilim]), label=ilim, lw=1.5, c=colors[k], source=True)

    except:
      print("\n\nLSST (Chang) p(z) not available.  Plotting whitebook p(z) estimate for %.4lf.\n\n" % ilims[ilim])

      pl.plot(zs, whitebook_pz(zs, ilims[ilim]), label=ilim, lw=1.5, c=colors[k])

  elg_pz       = elg_pz(interp=True)
  smooth_elgs  = convolve(elg_pz(zs), Box1DKernel(30))

  ##  Hildebrandt et al. (2009) -- measured dndz.
  midz, pz_H09 = getpz_H09()
  pz_H09       = interp1d(midz, pz_H09, kind='linear', bounds_error=False, fill_value=0.0)

  ##  BX:  'khaki'
  colors       = ['purple', 'lightskyblue', 'mediumseagreen',       'khaki',    'dodgerblue',     'darkgreen',             'r']
  labels       = ['DR12 QSOs', 'DESI ELGs',    'BM-dropouts', 'BX-dropouts', r'$u$-dropouts', r'$g$-dropouts', r'$r$-dropouts']

  zs, dVs      = dVols(zs, cosmo, params)

  zs           = np.concatenate([zs, np.array([zmax])])
  dVs          = np.concatenate([np.zeros(1),   dVs])

  ## BM / BX redshift 2 selection; https://arxiv.org/pdf/astro-ph/0401439.pdf 
  ## Hildebrandt et al. (2009) -- Gaussian approximation:  Gaussian(zs, 2.96, 0.24)
  curves       = [get_pz(zs, 'QSO'), smooth_elgs, Gaussian(zs, 1.70, 0.34), Gaussian(zs, 2.20, 0.32),\
                  pz_H09(zs), get_completeness(drop='g')(zs) * dVs, get_completeness(drop='r')(zs) * dVs]

  for pz, color, label in zip(curves, colors, labels):
    norm = np.sum(pz * dz)

    pl.plot(zs, pz / norm, color, label = label)

  pl.xlabel(r'$z$', fontsize=12)
  pl.ylabel(r'$p(z)$', fontsize=12)

  pl.xlim(0.0, 7.0)
  pl.ylim(0.0, 2.5)

  pl.legend(loc=1, ncol=3, frameon=False)

  plt.tight_layout()
  
  ##  pl.show()
  pl.savefig('plots/survey_pz.pdf')  
  
  print("\n\nDone.\n\n")
