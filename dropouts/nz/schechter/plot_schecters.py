import   numpy              as      np
import   pandas             as      pd
import   matplotlib         as      mpl
import   matplotlib.pyplot  as      plt

from     utils              import  pprint
from     nbar               import  comovdensity, projdensity
from     ilim               import  get_contamination

from     goldrush.specs     import  samplestats  as goldrush_stats
from     Malkan.specs       import  samplestats  as malkan_stats
## from  reddy              import  samplestats  as reddy_stats


plt.style.use('ggplot')
mpl.rc('text', usetex = True)

def plot_schechters():
    import  os
    import  pylab              as      pl
    import  matplotlib.pyplot  as      plt
    
    from    matplotlib.pyplot  import  text
    from    specs              import  samplestats
    from    selection_box      import  detection_bands
    from    utils              import  latexify


    latexify(fig_width=None, fig_height=None, columns=2, equal=False)
    
    colors = ['b', 'r', 'indigo']
    
    ##  Create figure.                                                                                                                        
    fig, axarray = plt.subplots(1, 1, sharey=False)
    fig.set_size_inches(6.5, 3.5)

    dm           = 0.25
    mlims        = np.arange(20.5, 28.0, dm)

    ##  The old fashioned way. 
    count        = 0

    samples      = [malkan_stats(), goldrush_stats(), goldrush_stats()]
    bands        = ['Malkan', 'g', 'r']
    labels       = [r'$u$-dropouts', r'$g$-dropouts', r'$r$-dropouts']
    
    stop         = 1

    for stats, band, label in zip(samples[:stop], bands[:stop], labels[:stop]):         
      zee        = stats[band]['z']   
      dzee       = 0.5

      alpha      = stats[band]['schechter']['alpha']
      Mstar      = stats[band]['schechter']['M_star']
      phi_star   = stats[band]['schechter']['phi_star']

      pnbars     = []

      for mlim in mlims:
        pnbar = projdensity(zee - dzee / 2., zee + dzee / 2., phi_star, Mstar, alpha, mlim=mlim, printit = True, completeness=None)
        pnbars.append(pnbar)
  
      pnbars  = np.array(pnbars)

      pl.semilogy(mlims, pnbars, '-', label=label, color=colors[count], alpha=0.5)

      ## Malkan u-dropouts contamination corrected. 
      if band == 'Malkan':
          pl.semilogy(mlims, (1. - get_contamination(mlims, zee=3, depth='W')) * pnbars, '^', label=label + ' (corrected)', color=colors[count],\
                      alpha=0.5, markersize=2)

      ## Gradient for magnification bias. 
      grads   = np.gradient(pnbars, dm)
      indexs  = (mlims / pnbars) * grads

      pl.semilogy(mlims, indexs, '--', c=str(colors[count]), dashes=[3,1])

      
      for i, n in enumerate(indexs):
          if i % 4 == 0:
            A   = pnbars[i] / mlims[i] ** n
          
            pl.semilogy(mlims, A * mlims ** n, '--', color='k', label='', alpha=0.5) 
      
      count     += 1

    pl.xlim(22.5,  27.0)
    pl.ylim(1.e0,  7.e4)

    pl.xlabel(r'$m_{\rm{UV}}$')
    pl.ylabel(r'$N(<m_{\rm{UV}})$ / deg$^2$')
      
    pl.legend(loc=2)

    pl.savefig('plots/schechters.pdf', bbox_inches='tight')


if __name__ == "__main__":
    import  collections


    print "\n\nWelcome to a Schechter fn. plotter.\n\n"

    depths = ['W']    ## ['W', 'D', 'UD'].              

    plot_schechters()

    ## pprint(stats)

    print "\n\nDone.\n\n"
