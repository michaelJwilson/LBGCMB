import  hmf
import  numpy              as      np
import  matplotlib         as      mpl
import  matplotlib.pyplot  as      plt
import  pylab              as      pl

from    hmf                import  cosmo
from    astropy.cosmology  import  FlatLambdaCDM
from    utils              import  latexify
from    params             import  get_params
from    scipy.interpolate  import  interp1d


latexify(columns=1, equal=True, fontsize=10, ggplot=True, usetex=True)

params =  get_params()
cosmo  = {'H0': 100. * params['h_100'], 'Om0': 0.3153, 'Ob0': 0.02242 / 0.6736 ** 2.}

mf     =  hmf.MassFunction(cosmo_params = cosmo)
models = [hmf.fitting_functions.ST, hmf.fitting_functions.Tinker08]

def get_stellar2halo(drop_band = 'u'):
  data           = np.loadtxt('../dat/stellar/smhm.txt')

  halo_mass      = data[:,0]               ## [Msun / h]

  dict           = {'u': data[:,1], 'g': data[:,2], 'r': data[:,3]}
  
  stellar_mass   =  dict[drop_band]
  stellar_mass   =  10. ** stellar_mass
  stellar_mass  *= (halo_mass / params['h_100'])  ## Mstellar [Msun]

  return  interp1d(stellar_mass, halo_mass, kind='linear', copy=True, bounds_error=True)

if __name__ == '__main__':
  stellar_masses   = 10. ** np.arange(7.5, 10.1, 0.1) 

  for drop_band, zee in zip(['u', 'g', 'r'], [3., 4., 5.]):
    stellar2halo   = get_stellar2halo(drop_band = drop_band)
    min_halomasses = stellar2halo(stellar_masses)

    result         = []

    for Mmin in np.log10(min_halomasses):
      mf.update(z  = zee, Mmin = Mmin, Mmax = np.log10(5.e14), hmf_model = models[0], delta_h = 200., n = 1.0)

      result.append(mf.ngtm[0])

    result = np.array(result)

    pl.loglog(stellar_masses, result, label=r'$z=$' + ' %.1lf' % zee)


ax     = pl.gca() 

ax.set_xlabel(r'$M_* \ [M_\odot]$')
ax.set_ylabel(r'$\bar n(> M_*) \ [(\rm{Mpc} / $$h$$)^{-3}]$')

pl.xlim(1.0e+7, 1.1e10)
pl.ylim(1.0e-4, 1.e00)

pl.legend(loc=3, ncol=1)

pl.savefig('../plots/stellarmass_nz.pdf', bbox_inches='tight')
