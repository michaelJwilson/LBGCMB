import  hmf
import  matplotlib         as      mpl
import  numpy              as      np
import  matplotlib.pyplot  as      plt

from    scipy.integrate    import  simps
from    scipy.interpolate  import  interp1d
from    hmf                import  cosmo
from    astropy.cosmology  import  FlatLambdaCDM

plt.style.use('ggplot')
mpl.rc('text', usetex = True)


if __name__ == "__main__":
    import  matplotlib         as      mpl
    import  matplotlib.pyplot  as      plt
    import  pylab              as      pl


    print("\n\nWelcome to massfn.\n\n")

    '''
    ## Load Martin's results.
    data = np.loadtxt('../dat/massfn_st_0.2000.dat')
    pl.loglog(data[:,0], data[:,1], label='Martin-ST')

    data = np.loadtxt('../dat/massfn_tin_0.2000.dat')
    pl.loglog(data[:,0], data[:,1], label='Martin-Tinker')
    '''
    ##
    ## Murray et al. hmf python code. 
    ## 
    ## ncosmo = FlatLambdaCDM(H0 = 73.0, Om0=0.25, Tcmb0 = 5.0, Ob0 = None)
    ## ncosmo = cosmo.Cosmology(cosmo_model   = ncosmo)

    for zee in [3.0, 4.0]: 
      pl.clf()

      ## Ok0': 0.
      for lcosmo, cosmo, lstyle in zip(['XXL', 'Planck (2018)'], [{'H0': 73.,   'Om0': 0.2500, 'Ob0': 0.05 * 0.25},\
                                                                  {'H0': 67.36, 'Om0': 0.3153, 'Ob0': 0.02242 / 0.6736 ** 2.}], ['-', '--']):
        mf = hmf.MassFunction(cosmo_params = cosmo)  

        for label, fitting, color in zip(['Sheth-Tormen', 'Tinker', 'Jenkins'], [hmf.fitting_functions.ST, hmf.fitting_functions.Tinker08, hmf.fitting_functions.Jenkins], ['b', 'r', 'g']):
          mf.update(z = 4.0, Mmin = 8., Mmax = np.log10(5.e14), hmf_model = fitting, delta_h = 200., n = 1.0)

          print  mf.parameter_values
          ## print  mf.parameter_info()

          pl.loglog(mf.m, mf.ngtm, color + lstyle, label=label + ': ' + lcosmo, alpha=0.6)

      ## Axis labels.
      ax   = pl.gca()

      ax.set_xlabel(r'$M_{\rm{min}} \ [M_\odot / h]$')
      ax.set_ylabel(r'$\bar n \ [(\rm{Mpc} / h)^{-3}]$')

      pl.xlim(1.e+11, 1.e14)
      pl.ylim(1.e-10, 1.e00)

      pl.legend(loc=3, ncol=2)

      pl.savefig('../plots/massfn_z%.0lf.pdf' % zee, bbox_inches='tight')
    
    print("\n\nDone.\n\n")
