import  matplotlib         as      mpl
import  matplotlib.pyplot  as      plt

from    params             import  get_params
from    utils              import  comoving_distance
from    zeldovich_Lmax     import  Lcutmax
from    pylab              import  rcParams
from    utils              import  latexify


latexify(fig_width=None, fig_height=None, columns=1, equal=True)

## rcParams['figure.figsize'] = 3.5, 3.5

params = get_params()

def angular2spatial(ell, redshift):
    K  = (ell + 1./2.)/comoving_distance(redshift)

    kz = 0.0
    k2 = K**2. + kz**2.

    return k2**0.5

def plot_angular2spatial():
  """
  Construct a plot showing the k = (L + 0.5) / chi(z) relation. 
  """
  import  matplotlib.pyplot  as plt
  import  numpy              as np
  import  pylab              as pl
  
  ## Plot Zel'dovich Lmax locus.                                                                                                              
  zs     = Lcutmax.keys()
  Lmaxes = np.array(Lcutmax.values())
   
  pl.plot(Lmaxes[:,0].astype('float32'), zs, 'k', lw=1.2)

  ells   = np.arange(500., 2510., 10.0)
  zs     = np.arange(0.50,   6.1,  0.1)

  levels = 0.05*np.arange(1, 25, 1)

  grid_ells, grid_zs  = np.meshgrid(ells, zs)

  grid_ks             = angular2spatial(grid_ells, grid_zs)
  CS                  = plt.contourf(grid_ells, grid_zs, grid_ks, levels=levels, cmap='tab20c')

  cbar                = plt.colorbar(CS)

  CS2                 = plt.contour(CS, colors='k', alpha=0.0)
  cbar.add_lines(CS2)

  plt.clabel(CS2, inline=True, fontsize=8, color='k', inline_spacing=10)

  cbar.ax.set_ylabel(r'$k \ [h \rm{Mpc}^{-1}]$')
  
  pl.xlim(500., 2000.)
  pl.ylim(0.9,    6.0)

  pl.xlabel(r"$L$")
  pl.ylabel(r"$z$")

  pl.savefig("plots/angular2spatial.pdf", bbox_inches='tight')


if __name__ == "__main__":
    plot_angular2spatial()
