import  matplotlib;                      matplotlib.use('PDF')

import  matplotlib.pyplot        as      plt
import  healpy                   as      hp
import  pylab                    as      pl
import  numpy                    as      np
import  ephem

from    mpl_toolkits.axes_grid1  import  make_axes_locatable


plt.style.use('ggplot')

def plot_mwd(ax, RA, Dec, org=0, title='', color='k', s=0.01, alpha=1.0, marker='.'):
    ''' 
    RA [0, 360), Dec [-90, 90] degrees are arrays of the same length.

    org is the origin of the plot, 0 or a multiple of 30 degrees in [0, 360).
    projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
    '''
    x       =  np.remainder(RA + 360. - org, 360.)  ##  Shift RA values
    ind     =  x > 180.
    x[ind] -=  360.                                 ##  Scale conversion to [-180, 180]
    x       = -x                                    ##  Reverse the scale: East to the left
    
    ax.scatter(np.radians(x), np.radians(Dec), color=color, s=s, rasterized=True, alpha=alpha, marker=marker)  

    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder(tick_labels + 360 + org, 360)

    ax.set_xticklabels(tick_labels, color='white')                        

    ax.set_title(title)
    ax.title.set_fontsize(15)

    ax.set_xlabel("Right Ascension [degs.]")
    ax.xaxis.label.set_fontsize(12)

    ax.set_ylabel("Declination [degs.]")
    ax.yaxis.label.set_fontsize(12)

    ax.grid(True)

def get_plane(input = ephem.Galactic, output = ephem.Equatorial):
    lon_array  = np.arange(0., 360., 1e-2)
    lat        = 0.0

    eq_array   = np.zeros((len(lon_array), 2))

    for i, lon in enumerate(lon_array):
        ga            = input(np.radians(lon), np.radians(lat))
        eq            = output(ga)

        eq_array[i]   = np.degrees(eq.get())

    RA  = eq_array[:,0]
    DEC = eq_array[:,1]

    return RA, DEC


if __name__ == '__main__':
    randoms    = hp.fitsfunc.mrdfits("/global/homes/m/mjwilson/desi/randoms/randoms.fits", hdu=1)

    ra         = randoms[0][::50]
    dec        = randoms[1][::50]

    fig        = plt.figure(figsize=(10, 5))
    ax         = fig.add_subplot(111, projection='mollweide', facecolor = 'gainsboro')


    plot_mwd(ax, ra, dec, org=90, title ='DESI footprint', color='dodgerblue', alpha=0.1)

    RA, DEC    = get_gplane()

    plot_mwd(RA, Dec, 180, color='k', alpha=1.0, s=0.5)
    
    pl.savefig("footprint.pdf")


