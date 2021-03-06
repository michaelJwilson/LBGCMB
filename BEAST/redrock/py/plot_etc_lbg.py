import  os
import  pickle
import  pylab                   as     pl
import  numpy                   as     np
import  matplotlib.pyplot       as     plt
import  matplotlib              as     mpl
import  matplotlib.colors       as     colors

from    astropy.table           import Table, vstack
from    utils                   import latexify
from    redrock.results         import read_zscan
from    mpl_toolkits.axes_grid1 import make_axes_locatable


dband      =     'r'
target     =   'lbg'
survey     =  'desi'
base_exp   =   1500.         ##  [Seconds].                                                                                                               
repeat     =       1
nexposures =      10         ##  Scaling,  not coadd.                                                                                                     

root       = os.environ['BEAST']
ddir       =      ''

##  Truth values.                                                                                                                                      
Truth     =  root + '/gal_maker/dat/Tables/galmaker-lbg-meta.txt'
Truth     =  Table.read(Truth, format='ascii', names=['OBJTYPE', 'SUBTYPE', 'REDSHIFT', 'Lya-EW', 'r'])
Truth     =  vstack([Truth] * repeat)

TRUEZS    =  np.array(Truth['REDSHIFT'].quantity).astype(np.float)
TRUEEWS   =  np.array(Truth['Lya-EW'].quantity).astype(str)
TRUEMAGS  =  np.array(Truth[dband].quantity).astype(np.float)

EWs       = list(set([str(x) for x in Truth['Lya-EW']]))

for EW in EWs:
    ##  Pickle it.                                                                                                                                                                   
    es             = pickle.load(open(root + '/redrock/pickle/%s/%s/exptimes.pkl' % (survey, target),     'rb'))
    es             = np.array(es).astype(np.float)
    es            /= 1500.  ##  Seconds to exposures.  

    ##  Grid of min. exposure time in hours for a given redshift, magnitude and equivalent width.                                                    
    ##  es = np.ma.masked_invalid(es)

    pl.clf()

    ##  And plot.                                                                                                                                         
    ##  latexify(fig_width=None, fig_height=None, columns=1, equal=True, ggplot=False, usetex=True)                                                      

    mpl.rcParams.update({'font.size': 14, 'font.family': 'monospace'})
    
    fig      = pl.gcf()

    cmap     = plt.cm.tab20c
    cmap     = colors.ListedColormap([cmap(i) for i in np.array([2, 1, 4, 5, 8, 9, 12, 13, 16, 17])[::1]])  ##  'viridis'                               

    ## Set bad to be white.                                                                                                                              
    cmap.set_bad(color = 'white', alpha = 1.)
    
    bounds   = np.arange(1, 10, 1)
    norm     = colors.BoundaryNorm(bounds, cmap.N, clip=False)

    ##  print(type(EW), np.array(Truth['Lya-EW']).astype(str), np.array(Truth['Lya-EW'])[np.array(Truth['Lya-EW']).astype(str) == EW])

    cut      = np.where((TRUEEWS == EW) & (es > 2.e4))
    pl.scatter(TRUEZS[cut], TRUEMAGS[cut], c='k', marker='x', alpha=0.2)

    print('Number with no redshifts in maximum exposure: %d' % len(es[cut]))

    cut      = np.where((TRUEEWS == EW) & (es < 2.e4))
    ##  plt.pcolormesh(xv, yv, zv, cmap=cmap, norm=norm) 
    pl.scatter(TRUEZS[cut], TRUEMAGS[cut], c=es[cut].astype(np.float), marker='^', cmap=cmap, norm=norm)
    
    plt.text(4.5, 26.25, EW+r'$\AA$')

    pl.xlim(np.unique(TRUEZS)[0], np.unique(TRUEZS)[-1])
    pl.ylim(np.unique(TRUEMAGS[np.isfinite(TRUEMAGS)])[0], np.unique(TRUEMAGS[np.isfinite(TRUEMAGS)])[-1])

    pl.xlabel(r'$z$')
    pl.ylabel(r'$r_{AB}$')

    ax      = pl.gca()
    divider = make_axes_locatable(ax)
    cax     = divider.append_axes("right", size="5%", pad=0.05)

    cb      = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds+0.5,\
                                        boundaries=bounds, format='%.2lf')

    cb.set_ticklabels(bounds, update_ticks=True)
    cb.set_label('Exposure time [1500 s]', rotation=270, labelpad=20)
        
    plt.tight_layout()

    pl.savefig(os.environ['BEAST'] + '/redrock/plots/%s/%s/expgrid_%s.pdf' % (survey, target, EW.replace('.', 'p')), bbox_inches='tight')
