import  numpy              as  np
import  pylab              as  pl
import  matplotlib.pyplot  as  plt

from    cosmo              import  cosmo
from    utils              import  latexify
from    schechter.gen_pz   import  peakz           as  _peakz
from    schechter.get_shot import  get_shot
from    schechter.get_pz   import  get_pz
from    nbar               import  comovdensity
from    get_bz             import  bz_callmodel
from    get_schechters     import  get_schechters
from    get_wschechters    import  get_wschechter
from    growth_rate        import  growth_rate
from    reddy.specs        import  samplestats     as reddy_stats
from    goldrush.specs     import  samplestats     as grush_stats
from    Malkan.specs       import  samplestats     as malkan_stats
from    reddy.pz           import  get_pz          as reddy_getpz
from    goldrush           import  completeness    as grush_completeness
from    Malkan             import  completeness    as malkan_completeness


latexify(columns=2, ratio=0.5, equal=False, fontsize=12, ggplot=True, usetex=True)

for color, band in zip(['y', 'b', 'g', 'r'], ['BX', 'u', 'g', 'r']):
    frac    =          True
    area    =         15000.

    deltav  =           400.    ##  [km / s].                                                                                                                    
    kmax    =           0.9

    pz      =  get_pz(band)
    peakz   =    _peakz(pz)

    fsky    =  area / 41253.
    sigp    =  (1. + peakz) * deltav  / cosmo.efunc(peakz) / 100.  ##  [Mpc / h].  
    
    dat     =  np.loadtxt('dat/rsd_%s_%.1lf_%.1lf_%.3lf.dat' % (band, kmax, fsky, deltav))

    print('%.3lf \t %.6lf' % (peakz, growth_rate(1. / (1. + peakz))))

    if frac:
      pl.plot(np.ones_like(dat[:,1][dat[:,0] < 0.5]) * peakz, 100. * dat[:,1][dat[:,0] < 0.5] / growth_rate(1. / (1. + peakz)),\
                           color + '_', alpha=0.5, markersize=20)
    
    else:
      pl.plot(np.ones_like(dat[:,1][dat[:,0] < 0.5]) * peakz, 100. * dat[:,1][dat[:,0] < 0.5], color + '_', alpha=0.5, markersize=20)

pl.xlim(2.0, 5.0)
pl.ylim(0.,  5.0)

pl.xlabel(r'$z$')
pl.ylabel(r'$(\sigma_f / f) \ [\%]$')
pl.legend(ncol=1, loc=2, frameon=False)

plt.tight_layout()

##  pl.show()
pl.savefig('plots/sigf.pdf')
