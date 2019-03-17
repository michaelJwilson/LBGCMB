import   numpy                   as      np
import   pandas                  as      pd
import   matplotlib              as      mpl
import   matplotlib.pyplot       as      plt
import   astropy.units           as      u

from     utils                   import  pprint
from     nbar                    import  comovdensity, projdensity
from     utils                   import  comoving_distance
from     cosmo                   import  cosmo
from     schechter.gen_pz        import  peakz            as  _peakz
from     schechter.get_shot      import  get_shot
from     schechter.get_pz        import  get_pz
from     goldrush                import  completeness     as  grush_completeness
from     goldrush.contamination  import  get_contamination
from     goldrush.ilim           import  get_nbarbymag
from     Malkan                  import  completeness     as  malkan_completeness
from     get_schechters          import  get_schechters
from     get_wschechters         import  get_wschechter
from     reddy.specs             import  samplestats      as  reddy_stats
from     goldrush.specs          import  samplestats      as  grush_stats
from     Malkan.specs            import  samplestats      as  malkan_stats
from     reddy.pz                import  get_pz           as  reddy_getpz


def plot_schechters(magbias=False, nointerloper=False):
    import  os
    import  pylab              as      pl
    import  matplotlib.pyplot  as      plt
    
    from    matplotlib.pyplot  import  text
    from    utils              import  latexify


    latexify(columns=2, equal=False, ggplot=True, fontsize=12, ratio=0.5)
        
    ##  Create figure.                                                                                                                        
    ##  fig, axarray = plt.subplots(1, 1, sharey=False)
    ##  fig.set_size_inches(6.5, 3.5)

    dm           = 0.25
    mlims        = np.arange(20.5, 28.0, dm)

    count        = 0
    pnbars       = []

    print('\n\nSolving for Quasars.')

    for mlim in mlims:
        pnbar    = projdensity(2., 3., None, None, None, mlim=mlim, type='qso')
        pnbars.append(pnbar)

    pnbars  = np.array(pnbars)
    pl.semilogy(mlims, pnbars, '-', label=r'$2 < z < 3$ QSO', color='c', alpha=0.5)
    
    ##  And save.
    np.savetxt(os.environ['LBGCMB'] + '/dropouts/schechter/dat/qso.dat', np.c_[mlims, pnbars], fmt='%.6lf', header='mlim   g/deg^2')

    ##  Now dropouts.
    samples      = [reddy_stats(), malkan_stats(), grush_stats(), grush_stats()]
    
    bands        = ['BX', 'Malkan', 'g', 'r']
    labels       = ['BX-dropouts', r'$u$-dropouts', r'$g$-dropouts', r'$r$-dropouts']
    colors       = ['y', 'b', 'g', 'indigo']

    stop         =  2

    ##  for stats, band, label, color in zip(samples[:stop], bands[:stop], labels[:stop], colors[:stop]):         
    for stats, band, label, color in zip(samples, bands, labels, colors):
      print('\n\nSolving for %s.' % label)
        
      dzee       = 0.7
      zee, alpha, Mstar, phi_star = get_schechters(stats, band)
      
      pnbars     = []
      
      for mlim in mlims:
        if band in ['g', 'r']:    
          ##  Apply completeness correction to Schechter estimate.   
          CC     = grush_completeness.get_completeness(band)  

          pnbar  = projdensity(zee - dzee / 2., zee + dzee / 2., phi_star, Mstar, alpha, mlim=mlim, type='app', printit=True,\
                               completeness= lambda zz: CC(zz))
        else:
         pnbar   = projdensity(zee - dzee / 2., zee + dzee / 2., phi_star, Mstar, alpha, mlim=mlim, type='app', printit=True)   

        pnbars.append(pnbar)
      
      pnbars = np.array(pnbars)

      pl.semilogy(mlims, pnbars, '-', label=label, color=color, alpha=0.5)
            
      ##  And save.                                                                                                                                 
      np.savetxt(os.environ['LBGCMB'] + '/dropouts/schechter/dat/%sDrop.dat' % band, np.c_[mlims, pnbars], fmt='%.6lf', header='mlim   g/deg^2')

      if band == 'BX':                                                                                                                              
        for mlim in np.arange(22.5, 26.0, 0.5):                                                                                                                                                                                           
          stats = reddy_stats(mlim)                                                                                                          

          pl.semilogy(mlim, stats['BX']['nbar'],       'y^', markersize=3)
          pl.semilogy(mlim, stats['BX']['nbar_noint'], 'ys', markersize=3)
          
      if band == 'Malkan':
        ##  Only 0.5 sampling available. 
        nbars  = []
        crates = []

        for mlim in mlims[::2]:
          stats = malkan_stats(mlim)
          nbar  = stats['Malkan']['nbar'] 

          nbars.append(nbar)
          crates.append(get_contamination(mlim, round(stats['Malkan']['z']), depth='D'))

        nbars = np.array(nbars)
        pl.semilogy(mlims[::2], nbars, 'b^', markersize=3)
          
        crates = np.array(crates)
        pl.semilogy(mlims[::2], nbars * (1. - crates), color=color, marker='s', markersize=3, lw=0.)
        
      if band == 'g':
        magbins, pnbar, counts = get_nbarbymag(band, depth='W', printit=False)
        pl.semilogy(magbins[:-1:50], pnbar[::50], color=color, marker='^', markersize=3, lw=0.)

        rate = get_contamination(magbins[:-1:50], round(stats[band]['z']), depth='W')
        pl.semilogy(magbins[:-1:50], pnbar[::50] * (1. - rate), color=color, marker='s', markersize=3, lw=0.)

      if band == 'r':
        magbins, pnbar, counts = get_nbarbymag(band, depth='D', printit=False)
        pl.semilogy(magbins[:-1:50], pnbar[::50], color=color, marker='^', markersize=3, lw=0.)

        rate = get_contamination(magbins[:-1:50], round(stats[band]['z']), depth='D')
        pl.semilogy(magbins[:-1:50], pnbar[::50] * (1. - rate), color=color, marker='s', markersize=3, lw=0.)
      
      if magbias:
        ##  Gradient for magnification bias. 
        grads   = np.gradient(pnbars, dm)
        indexs  = (mlims / pnbars) * grads

        pl.semilogy(mlims, indexs, '--', c=str(color), dashes=[3,1])

        for i, n in enumerate(indexs):
          if i % 4 == 0:
            A = pnbars[i] / mlims[i] ** n
          
            pl.semilogy(mlims, A * mlims ** n, '--', color='k', label='', alpha=0.5) 
       
    pl.xlim(22.4,  27.1)
    pl.ylim(1.e0,  7.e4)

    pl.xlabel(r'$m$')
    pl.ylabel(r'$N(<m)$ / deg$^2$')

    pl.legend(loc=4, ncol=2)
    
    '''
    ##  Line Schecter functions. 
    ax  = pl.gca()
    axx = ax.twiny()

    ##  axx.set_xlabel(r'log$_{10}$\{$L_{\rm{min}} / $(ergs$/s$)\}')

    ##  Lyman-alpha. 
    ##  Lower limits on LF integral.
    logLmins = np.arange(40.0, 44.0, 0.1)

    pnbars   = []

    for logLmin in logLmins:                                                                                 
      pnbar  = projdensity(2.2, 3.5, None, None, None, mlim=logLmin, type='lya', printit=True)                                                                    
      pnbars.append(pnbar)                                                                                                                       

    pnbars   = np.array(pnbars)                                                                                                                        
    axx.semilogy(logLmins, pnbars, '-', label=r'$2.2 < z < 3.5$ Ly-$\alpha$', color='m', alpha=0.5, zorder=-1)
    '''
    '''
    ##  ... and Comparat OII.
    print('\n\nSolving for Comparat OII.\n\n')

    logSmins  = np.arange(-21., -14.0, 0.1)[::-1]
    pnbars    = []

    for logSmin in logSmins:
      pnbar   = projdensity(0.6, 1.6, None, None, None, mlim=logSmin, type='oii', printit=False, app_linelim=True)
      pnbars.append(pnbar)

    pnbars    = np.array(pnbars)

    axx.semilogy(logSmins, pnbars, '-', label=r'$0.6 < z < 1.6$ OII', color='c', alpha=0.5, zorder=-1)

    ##  axx.set_xlim(44., 39.)
    
    pl.legend(loc=2)
    '''
    
    ##  pl.show()
    pl.savefig('plots/schechters.pdf', bbox_inches='tight')


if __name__ == "__main__":
    import  collections


    print("\n\nWelcome to a Schechter fn. plotter.\n\n")

    plot_schechters()

    ## pprint(stats)

    print("\n\nDone.\n\n")
