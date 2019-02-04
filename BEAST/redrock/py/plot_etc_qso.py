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
target     =   'qso'
survey     =  'desi'
base_exp   =   1500.         ##  [Seconds].                                                                                                               
repeat     =       1
nexposures =      10         ##  Scaling,  not coadd.                                                                                                     

root       = os.environ['BEAST']

##  Truth values.                                                                                                                                            
names     = ['TARGETID', 'OBJTYPE', 'SUBTYPE', 'TEMPLATEID', 'SEED', 'REDSHIFT', 'MAG', 'MAGFILTER', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2']

Truth     = os.environ['BEAST'] + '/gal_maker/dat/Tables/galmaker-qso-meta.txt'
Truth     = Table.read(Truth, format='ascii', names=names)

Truth     = Table([Truth['TARGETID'], Truth['OBJTYPE'], Truth['SUBTYPE'], Truth['REDSHIFT'], Truth['MAG'], Truth['MAGFILTER']],\
                   names=('TARGETID', 'OBJTYPE', 'SUBTYPE', 'REDSHIFT', 'MAG', 'MAGFILTER'))

Truth     = vstack([Truth] * repeat)
  
TRUEZS         =  np.array(Truth['REDSHIFT'].quantity).astype(np.float)
TRUEEWS        =  np.zeros_like(Truth['REDSHIFT'].quantity).astype('str')  ##  np.array(Truth['EWOII'].quantity).astype(str)                                 
TRUEMAGS       =  np.array(Truth['MAG'].quantity).astype(np.float)

EWs            = [str(x) for x in ['0.0']]
results        =  dict(zip(EWs, [[] for x in EWs]))

for EW in EWs:
    ##  Pickle it.                                                                                                                                     
    es = pickle.load(open(root + '/redrock/pickle/%s/%s/exptimes.pkl' % (survey, target), 'rb'))
    es /= 1500.  ##  S to exposures. 

    ##  Grid of min. exposure time in hours for a given redshift, magnitude and equivalent width.                                                    
    es = np.ma.masked_invalid(es)

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

    cut      = np.where((TRUEEWS == EW) & (es > 2.e4))
    pl.plot(TRUEZS[cut], TRUEMAGS[cut], 'kx', alpha=0.2)
    
    print('Number with no redshifts in maximum exposure: %d' % len(es[cut]))

    cut      = np.where((TRUEEWS == EW) & (es < 2.e4))
    ## plt.pcolormesh(Truth['REDSHIFT'].quantity[cut], Truth[dband].quantity[cut], color=es[cut], cmap=cmap, norm=norm)
    pl.scatter(TRUEZS[cut], TRUEMAGS[cut], c=es[cut].astype(np.float), marker='^', cmap=cmap, norm=norm)

    plt.text(3.00, 22.0, EW + r'$\AA$')

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

    pl.savefig(os.environ['BEAST'] + '/redrock/plots/%s/%s/expgrid_%s.pdf' % (survey, target, str(EW).replace('.', 'p')), bbox_inches='tight')
