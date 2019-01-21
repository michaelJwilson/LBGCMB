import  matplotlib         as      mpl
import  numpy              as      np
import  glob
import  os
import  pylab              as      pl
import  matplotlib.pyplot  as      plt

from    utils              import  latexify


def plot_wt():
    '''
    Makes a 2 panel plot of w(theta) and b(r).
    '''
    tlo, thi = 30.0, 1000.0
    arcsec   = np.pi / 180. / 3600.
    
    latexify(fig_width=None, fig_height=None, columns=1, equal=True, fontsize=10)
    
    ##  First plot w(theta) for z=3 and 4.    
    cc    = 'darkblue'

    ## fname = os.environ['LBGCMB'] + '/Hildebrandt09/wtheta/H09_udrop_r24.5.txt'
    fname = os.environ['LBGCMB'] + '/Hildebrandt09/wtheta/wtheta_udropouts_r23p0t24p5_Masim.dat'

    wt    = np.loadtxt(fname)
    ww    = np.nonzero(wt[:,0] * 60. > tlo)[0]

    pl.errorbar(wt[ww,0] * 60., wt[ww,1], yerr=wt[ww,2], fmt='s', color=cc, label=r'$u$-drops (H09)', markersize=3)

    ##  Model
    wt = np.loadtxt("../dat/summary_v2/drop_0.2500_12.33_0.40.wt")
    pl.plot(wt[:,0] * 60., wt[:,1], c=cc, label='$z \simeq 3$ HOD', alpha=0.5)

    cc    = 'g'

    ##  Goldrush correlation fn. measurements.  To be found.                                                                  
    fname = os.environ['LBGCMB'] + '/dropouts/goldrush/cats/acfs/acf_converted_acf_W_combined_gri_23.50_24.50.data'
    wt    = np.loadtxt(fname)
    ww    = np.nonzero(wt[:,0] > tlo)[0]

    pl.errorbar(wt[ww,0], wt[ww,1], yerr=wt[ww,2], fmt='s', color=cc, label=r'$g$-drops (O17)', markersize=3)
    
    ##  Model
    wt = np.loadtxt("../dat/summary_v2/drop_0.2000_12.60_0.60.wt")
    pl.plot(wt[:,0]*60., wt[:,1], '-', color=cc, label='$z \simeq 4$ HOD', alpha=0.5)
    
    pl.xlim(20., thi)
    pl.ylim(1e-2,0.5)

    pl.xscale('log')
    pl.yscale('log')

    pl.legend(framealpha=0.2, loc=3, ncol=1)

    pl.xlabel(r'$\theta$  [arcsec]')
    pl.ylabel(r'$w(\theta)$')
    
    plt.savefig("../plots/wt.pdf", bbox_inches='tight')

def plot_bias():
    ##  Now plot the bias.
    latexify(fig_width=None, fig_height=None, columns=1, equal=True, fontsize=10)
    
    ##  z = 3.  
    br = np.loadtxt('../dat/summary_v2/drop_0.2500_12.33_0.40.brr')
    pl.plot(br[:,0], br[:,1], 's', c='darkblue', label=r'$b_{gg}(z \simeq 3)$', markersize=4)
    pl.plot(br[:,0], br[:,1], '-', c='darkblue', alpha=0.4)
    
    pl.plot(br[:,0], br[:,2], 'd', c='darkblue', label=r'$b_{gm}(z \simeq 3)$', markersize=4)
    pl.plot(br[:,0], br[:,2], '--', c='darkblue', dashes=[3,1], alpha=0.4)

    ## z = 4.
    br = np.loadtxt('../dat/summary_v2/drop_0.2000_12.60_0.60.brr')
    pl.plot(br[:,0], br[:,1], 's', color='darkgreen', label=r'$b_{gg}(z \simeq 4)$', markersize=4)
    pl.plot(br[:,0], br[:,1], '-', color='darkgreen', alpha=0.4)

    pl.plot(br[:,0], br[:,2], 'd',  color='darkgreen', label=r'$b_{gm}(z \simeq 4)$', markersize=4)
    pl.plot(br[:,0], br[:,2], '--', color='darkgreen', dashes=[3,1], alpha=0.4)

    pl.legend(framealpha=0.2)
    
    pl.xlabel(r'$r \ \ [h^{-1}{\rm Mpc}]$')
    pl.ylabel(r'$b(r)$')

    pl.xscale('log')
    pl.yscale('linear')

    pl.xlim(1.0, 40.)
    pl.ylim(3.7, 9.7)
    
    plt.savefig("../plots/br.pdf", bbox_inches='tight')


if __name__=="__main__":
  print('\n\nWelcome to a (plotter) of the angular correlation fn.\n\n')

  plot_wt()
    
  ## plot_bias()

  print('\n\nDone.\n\n')
