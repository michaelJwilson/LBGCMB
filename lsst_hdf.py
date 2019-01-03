import  numpy   as      np
import  pandas  as      pd
import  pylab   as      pl

from    utils   import  latexify


latexify(fig_width=None, fig_height=None, columns=2, equal=False)

"""
Compare the LSST N(i_AB) form assumed to measurements of the Hubble (Ultra) Deep Field (Metcalfe et al.) and Capak ++.
"""

def lsst_Ng(ilim):
    result = 46.*10.**(0.3*(ilim - 25.))  ## LSST science book (3.7); 46 galaxies per arcmin^2 for ilim = 25.

    return result * 60.**2.               ## galaxies per sq. deg.

def plot_Metcalfe():
    data  = pd.read_table("dat/lsst_nz/Metcalfe00.dat", skiprows=1, sep=r"\s*", names=['ilo', 'ihi', 'NorNraw', 'NorNgal', 'SouNraw', 'SouNgal'], engine='python')
    
    print data

    data['ilo'] += 0.37      ## Vega to AB conversion (http://www.astronomy.ohio-state.edu/~martini/usefuldata.html)                   
    data['ihi'] += 0.37      ## Vega to AB conversion 

    pl.semilogy(data['ilo'], 2.*data['NorNgal'], '*', markersize=3, label='HDF-N Metcalfe \'00')
    pl.semilogy(data['ilo'], 2.*data['SouNgal'], '*', markersize=3, label='HDF-S Metcalfe \'00')

    ## newdata = [[data['ilo'][2*i], data['NorNgal'][2*i] + data['NorNgal'][2*i+1]] for i, x in enumerate(data['ilo'][::2])]
    ## newdata = np.array(newdata)
    
    ## pl.semilogy(newdata[:,0], newdata[:,1], '*', markersize=3, label='HDF-N Metcalfe \'00')

    ## LSST
    los     = np.arange(18., 29., 1)
    pl.semilogy(los, lsst_Ng(1. + los) - lsst_Ng(los), label='LSST')

def plot_Nvilim():
    data = pd.read_table("dat/lsst_nz/capak07_Tab14.dat", skiprows=1, sep=r"\s*", names=['ilim', 'N', 'Poisson'], engine='python')
    pl.semilogy(data['ilim'], data['N'], '*', markersize=3, label='Capak \'07')

    plot_Metcalfe()

    pl.xlabel(r'$i_{AB}$')
    pl.ylabel(r'Number / deg$^2$ / mag')

    pl.xlim(22.,     28.)
    pl.ylim(10.**4., 10.**6.)

    pl.yscale('log')

    pl.legend(loc=2)

    pl.savefig('plots/lsst_nz_HDF.pdf', bbox_inches='tight')


if __name__ == '__main__':
    plot_Nvilim()
