import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    utils              import  latexify
from    scipy.interpolate  import  interp1d
from    whitebook_pz       import  Ng, lsst_chang_pz
from    schmittfull_bz     import  boss_lrg, desi_bgs, desi_lrg, desi_elg, desi_qso
from    growth_rate        import  growth_factor


if __name__ == '__main__':
    print('\n\nWelcome to DESI p(z).\n\n')

    ##  Load available BOSS / DESI samples. 
    bgs        = np.loadtxt('dat/nz_BGS.dat')
    lrgs       = np.loadtxt('dat/nz_LRG.dat')
    elgs       = np.loadtxt('dat/nz_ELG.dat')
    qsos       = np.loadtxt('dat/nz_QSO.dat')

    boss       = np.loadtxt('dat/nz_LRG.dat')
    boss[:,2]  = boss[:,3]
    
    zs         = np.arange(0.0, 5.0, 0.01)

    total_ns   = np.zeros_like(zs)
    total_cs   = np.zeros_like(zs)

    mean_bias  = np.zeros_like(zs)

    colors     = []

    latexify(columns=2, equal=False, fontsize=10, ratio=None, ggplot=False)

    samples    = [boss,     bgs,      lrgs,    elgs,      qsos]
    labels     = ['BOSS',   'BGS',   'LRGs',   'ELGs',   'QSOs']
    biases     = [boss_lrg, desi_bgs, desi_lrg, desi_elg, desi_qso]
    zmaxes     = [1.0, 0.5, 1.3, 1.8, 5.0]
    
    for (sample, label, bz, zmax) in zip(samples, labels, biases, zmaxes):
        csum       = np.cumsum(sample[:,2])
        icsum      = interp1d(np.concatenate([np.zeros(1), sample[:,1]]), np.concatenate([np.zeros(1), csum]), kind='linear',\
                              fill_value=(csum[0], csum[-1]), bounds_error=False)

        isum       = interp1d(np.concatenate([np.zeros(1), sample[:,1]]), np.concatenate([np.zeros(1), sample[:,2]]), kind='linear',\
                              fill_value=(0.0, 0.0), bounds_error=False)

        ns         =  isum(zs)
        cs         = icsum(zs)

        total_ns  += ns
        total_cs  += cs

        interim    = ns * bz(zs, force=True)
        interim[zs > zmax] = 0.0  

        mean_bias += interim
        
        p          = pl.plot(zs, cs, label=label)
        colors.append(p[0].get_color())

    ##  -- LSST --  
    depth        = 25.3                          ##  Y10  

    LSSTY10_nbar = Ng(depth, fmask=0.12)         ##  Per sq. deg.
    LSSTY10_pz   = lsst_chang_pz(zs, ilim=depth, source = False)

    LSSTY10_nz   = LSSTY10_nbar * LSSTY10_pz
    LSSTY10_cz   = np.cumsum(LSSTY10_nz)

    print('\n\nn_eff = %.4lf [per sq. deg.].' % LSSTY10_nz[-1])

    p = pl.plot(zs, LSSTY10_cz / 1.e4, 'magenta', label='LSST Y10 [$10^4$]')
    colors.append(p[0].get_color())

    ## Prep. mean bias.                                                                                                                                  
    mean_bias     /= total_ns
    mean_bias      = interp1d(zs, mean_bias, kind='linear', fill_value=(mean_bias[0], mean_bias[-1]), bounds_error=False)

    p              = pl.plot(zs,  total_cs, 'k-',  label='LOWZ + CMASS + DESI')
    colors.append(p[0].get_color())

    pl.xlim(3.0,    0.0)
    pl.ylim(0.0, 3.50e3)

    pl.xlabel(r'$z_{\rm{max}}$')
    pl.ylabel(r'$N(< z_{\rm{max}}) \ \ [\rm{gals. \ / \ deg^2}]$')

    ##  Second b(z) axis.                                                                                   
    ax  = pl.gca()
    axx = ax.twinx()

    ##  LSST sources
    ##  biases.append(lambda z: 1. + z)

    ##  LSST lenses
    biases.append(lambda z:  0.95 / growth_factor(1. / (1. + z)))

    biases.append(mean_bias)

    zmaxes.append(3.5) ## LSST bias limit
    zmaxes.append(3.5) ## DESI mean bias limit

    for i, (bz, zmax) in enumerate(zip(biases, zmaxes)):
        axx.plot(zs[zs < zmax], bz(zs[zs < zmax]), c=colors[i], alpha=0.45, dashes=[3,1], label='')

    axx.set_ylabel(r'$b(z)$')
    axx.set_ylim(0.0, 4.0)

    lines, labels = ax.get_legend_handles_labels()

    lines.insert(2, plt.Line2D([],[], alpha=0))
    labels.insert(2,'')

    lines.insert(5, plt.Line2D([],[], alpha=0))
    labels.insert(5,'')

    leg = ax.legend(lines, labels, loc=1, ncol=3, fancybox=True)
    leg.get_frame().set_alpha(1.00)

    pl.savefig('plots/desi_pz.pdf')

    ##  Save zs, cumulative N for all DESI (+ BOSS) samples, and effective bias.  
    np.savetxt('dat/effective_desi.dat', np.c_[zs[1:], total_cs[1:], mean_bias(zs[1:])], fmt='%.6lf \t')
    
    print('\n\nDone.\n\n')
