import  numpy              as  np
import  pylab              as  pl
import  pandas             as  pd
import  matplotlib.pyplot  as  plt


def plot_emlines(ymax=None):
    df     = pd.read_csv('./dat/sdss7_emlines.dat', names=['wave', 'gal', 'qso', 'line'], header=3, sep='\s+')
    df     = df.dropna()  

    ## Select according to list
    df     = df.loc[df['line'].isin(['Lya', 'OII', 'Ha', 'HB', 'CIV', 'Mg', 'Na', 'MgII'])]

    for col in ['gal', 'qso', 'wave']:
        df[col] = df[col].astype(np.float32)
        
    gf       = df

    lines    = gf['wave'].values

    Txlabels = {'Lya': r'Ly-$\alpha$', 'OII': r'OII', 'Ha': r'H$\alpha$', 'HB': r'H$\beta$', 'CIV': 'CIV', 'Mg': 'Mg', 'Na': 'Na', 'MgII': 'MgII'}
    labels   = gf['line'].values

    ratios   = gf['gal'].values
    ratios  /= ratios.max()

    plotted  = []
    shift    = 1.0

    if ymax is None:
        axes     = pl.gca()
        ylims    = axes.get_ylim()
        ymax     = ylims[1]

    props        = dict(boxstyle='round', facecolor='gainsboro', alpha=0.8, lw=0.0)
    
    for i, line in enumerate(lines):
        pl.axvline(x=line, ymin=0.0, ymax=1.0, alpha=0.7, c='k', lw=0.4)
  
        plt.text(line, ymax * (0.86 + 0.05 * shift), Txlabels[labels[i]],\
                 rotation=90, horizontalalignment='center', verticalalignment='center', bbox=props)

        plotted.append(labels[i])
  
        shift *= -1.0


if __name__ == "__main__":
    plot_emlines()

    elg      = np.loadtxt("./dat/elg-avg-line-ratio.dat")
    pl.plot(elg[:,0], elg[:,1], 'k')

    xs       = np.arange(1e3, 1e4,  10.)
    pl.plot(xs, np.zeros_like(xs), 'k')

    pl.xlim(1e3, 7.e3)

    pl.xlabel(r'$\AA$')
    pl.ylabel(r'$F_{\nu}$ [M Jy]')
    
    pl.savefig("plots/elg.pdf")
