import   numpy              as      np
import   pandas             as      pd
import   matplotlib         as      mpl
import   matplotlib.pyplot  as      plt

from     utils              import  pprint
from     nbar               import  comovdensity, projdensity
from     ilim               import  get_contamination

from     goldrush.specs     import  samplestats  as goldrush_stats
from     Malkan.specs       import  samplestats  as malkan_stats
from     reddy              import  samplestats  as reddy_stats


def plot_schechters(magbias=False, nointerloper=False):
    import  os
    import  pylab              as      pl
    import  matplotlib.pyplot  as      plt
    
    from    matplotlib.pyplot  import  text
    from    specs              import  samplestats
    from    selection_box      import  detection_bands
    from    utils              import  latexify


    latexify(fig_width=None, fig_height=None, columns=2, equal=False)
        
    ##  Create figure.                                                                                                                        
    fig, axarray = plt.subplots(1, 1, sharey=False)
    fig.set_size_inches(6.5, 3.5)

    dm           = 0.25
    mlims        = np.arange(20.5, 28.0, dm)

    count        = 0
    pnbars       = []

    ##  QSOs first. 
    for mlim in mlims:
        pnbar    = projdensity(2., 3., None, None, None, mlim=mlim, type='qso')
        pnbars.append(pnbar)

    pnbars  = np.array(pnbars)
    pl.semilogy(mlims, pnbars, '-', label=r'$2 < z < 3$ QSO', color='c', alpha=0.5)

    ##  Now dropouts.
    samples      = [reddy_stats(), malkan_stats(), goldrush_stats(), goldrush_stats()]
    bands        = ['BX', 'Malkan', 'g', 'r']
    labels       = ['BX-dropouts', r'$u$-dropouts', r'$g$-dropouts', r'$r$-dropouts']
    colors       = ['y', 'b', 'r', 'indigo']

    stop         =  2

    ## for stats, band, label, color in zip(samples[:stop], bands[:stop], labels[:stop], colors[:stop]):         
    for stats, band, label, color in zip(samples, bands, labels, colors):
      zee        = stats[band]['z']   
      dzee       = 0.5

      alpha      = stats[band]['schechter']['alpha']
      Mstar      = stats[band]['schechter']['M_star']
      phi_star   = stats[band]['schechter']['phi_star']

      pnbars     = []

      for mlim in mlims:
        pnbar = projdensity(zee - dzee / 2., zee + dzee / 2., phi_star, Mstar, alpha, mlim=mlim, type='app', printit = True, completeness=None)
        pnbars.append(pnbar)
  
      pnbars  = np.array(pnbars)

      pl.semilogy(mlims, pnbars, '-', label=label, color=color, alpha=0.5)

      ##  Malkan u-dropouts contamination corrected. 
      if (band == 'Malkan') & (nointerloper):
          pl.semilogy(mlims, (1. - get_contamination(mlims, zee=3, depth='W')) * pnbars, '^', label=label + ' (corrected)', color=color,\
                      alpha=0.5, markersize=2)

      if band == 'Malkan':
          ##  Only 0.5 sampling available. 
          for mlim in mlims[::2]:
              stats = malkan_stats(mlim)
              nbar  = stats['Malkan']['nbar'] 

              pl.plot(mlim, nbar, 'b^', markersize=3)

      if band == 'BX':
          for mlim in np.arange(22.5, 26.0, 0.5):
              stats = reddy_stats(mlim)
              nbar  = stats['BX']['nbar']

              pl.plot(mlim, nbar, 'y^', markersize=3)

      if magbias:
        ##  Gradient for magnification bias. 
        grads   = np.gradient(pnbars, dm)
        indexs  = (mlims / pnbars) * grads

        pl.semilogy(mlims, indexs, '--', c=str(color), dashes=[3,1])

        for i, n in enumerate(indexs):
          if i % 4 == 0:
            A = pnbars[i] / mlims[i] ** n
          
            pl.semilogy(mlims, A * mlims ** n, '--', color='k', label='', alpha=0.5) 

    pl.xlim(22.5,  27.0)
    pl.ylim(1.e0,  7.e4)

    pl.xlabel(r'$m_{5\sigma}$')
    pl.ylabel(r'$N(<m_{5\sigma})$ / deg$^2$')
      
    pl.legend(loc=4, ncol=2)

    pl.savefig('plots/schechters.pdf', bbox_inches='tight')


if __name__ == "__main__":
    import  collections


    print "\n\nWelcome to a Schechter fn. plotter.\n\n"

    depths = ['W']    ## ['W', 'D', 'UD'].              

    plot_schechters()

    ## pprint(stats)

    print "\n\nDone.\n\n"
