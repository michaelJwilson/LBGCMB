import matplotlib
import  os
import  numpy              as      np
import  pylab              as      pl
import  matplotlib.pyplot  as      plt


plt.style.use('ggplot')

if __name__ == '__main__':
    colors = ['b', 'r', 'indigo']

    ## Effective area
    area   =  876        ## [arcmin^2]
    area  /=  60. * 60.  ## [deg^2]


    fig, axarray = plt.subplots(1, 1, sharey=False)
    '''
    ## Observed
    data    = np.loadtxt('dat/Table2.dat')
    counts  = np.cumsum(data[:,2])
    counts /= area

    axarray.semilogy(data[:,1], counts, 'b-')
    '''
    ## Schechter estimate 
    root   = os.environ['LBGCMB']

    ## Malkan u-drops.
    data   = np.loadtxt(root + "/dropouts/nz/schechter/dat/schechter_estimate_Malkan_dropouts.txt")            
    axarray.semilogy(data[:,0], data[:,1], 'b', label='u-dropouts')

    Total  = np.copy(data[:,1])

    ## Goldrush g-dropouts.
    data   = np.loadtxt(root + "/dropouts/nz/schechter/dat/schechter_estimate_g_dropouts.txt")
    axarray.semilogy(data[:,0], data[:,1], 'g', label='g-dropouts')

    Total += data[:,1]
            
    ## Goldrush r-dropouts. 
    data   = np.loadtxt(root + "/dropouts/nz/schechter/dat/schechter_estimate_r_dropouts.txt")
    axarray.semilogy(data[:,0], data[:,1], 'r', label='r-dropouts')

    Total += data[:,1]
    axarray.semilogy(data[:,0], Total, 'k', label='Total')

    axarray.set_xlim(22.5,       25.5)                                                                                                         
    axarray.set_ylim([0.1, 5.e3])

    axarray.set_xlabel(r'5$\sigma$ mag. limit, $m$')
    axarray.set_ylabel(r'$N_{\rm{g}}(<m)$ / deg$^2$')

    axarray.legend(loc=2)

    pl.savefig('plots/ilims.pdf')
    
    print('\n\nDone.\n\n')
