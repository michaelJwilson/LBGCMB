import  ephem 
import  numpy                    as      np
import  pylab                    as      pl
import  healpy                   as      hp
import  pylab                    as      pl
import  matplotlib.pyplot        as      plt
import  matplotlib.animation     as      animation

from    footprint                import  plot_mwd, get_plane
from    utils                    import  latexify
from    mpl_toolkits.axes_grid1  import  make_axes_locatable
from    overlap                  import  prep_masks, plots


plt.style.use('ggplot')

latexify(columns=2, equal=False, fontsize=10, ratio=None, ggplot=True)

def hours2degs(input, dec=False):
    input   = input.split(':')

    ## One hr is 15 deg., One min is 1/60 of an hour. 
    scale   = [15., 15. / 60., 1.]
    
    if dec:
      scale = [1., 1./60., 1./60./60.]

    out     = np.array([np.float(x) * scale[i] for i, x in enumerate(input)])
    out     = np.sum(out)

    return  out

def cut(imap, pp, tt, lora, hira, declo, dechi):
    new               = np.ones_like(imap)
    
    new[pp < lora]    = 0.0
    new[pp > hira]    = 0.0
    
    new[tt < declo]   = 0.0
    new[tt > dechi]   = 0.0    

    return  new


if __name__ == '__main__':
    print('\n\nWelcome.\n\n')

    nside   =  256
    npix    =  hp.pixelfunc.nside2npix(nside)
    
    imap    =  np.arange(npix)

    t,  p   =  hp.pixelfunc.pix2ang(nside, imap)
    tt, pp  =  hp.Rotator(coord=['C', 'C'])(t, p)

    ##  dec-like
    tt      =  np.pi / 2. - tt
    tt     *=  180. / np.pi

    ##  ra-like
    pp     +=  np.pi 
    pp     *=  180. / np.pi

    print  tt.min(), tt.max()
    print  pp.min(), pp.max()
    
    ##  08:30 <= RA <= 15:00  &  -2 <= dec <= 5    (680 deg2)
    one     =  cut(imap, pp, tt, hours2degs('08:30:00'), hours2degs('15:00:00'), -2.,   5.)

    print hours2degs('08:30:00'), hours2degs('15:00:00')

    ##  22:00 <= RA <= 02:40  &  -1 <= dec <= 7    (630 deg2)  
    two     =  cut(imap, pp, tt, hours2degs('22:00:00'), hours2degs('24:00:00'), -1.,   7.)
    twotwo  =  cut(imap, pp, tt, hours2degs('00:00:00'), hours2degs('02:40:00'), -1.,   7.)

    print hours2degs('22:00:00'), hours2degs('24:00:00')
    print hours2degs('00:00:00'), hours2degs('02:40:00')

    ##  13:20 <= RA <= 16:40  &  42 <= dec <= 44.5 ( 90 deg2) 
    three   =  cut(imap, pp, tt, hours2degs('13:20:00'), hours2degs('16:40:00'), 42., 44.5)

    print hours2degs('13:20:00'), hours2degs('16:40:00')

    ## Take the union.
    all     =  one + two + twotwo + three

    ##  Sanity check on patch overlap
    print np.unique(all)

    print
    print
    print  np.sum(one)      *  hp.pixelfunc.nside2pixarea(nside, degrees = True)
    print  np.sum(two)      *  hp.pixelfunc.nside2pixarea(nside, degrees = True)
    print  np.sum(twotwo)   *  hp.pixelfunc.nside2pixarea(nside, degrees = True)
    print  np.sum(three)    *  hp.pixelfunc.nside2pixarea(nside, degrees = True)
    print
    print 
    print  np.sum(all)      *  hp.pixelfunc.nside2pixarea(nside, degrees = True)

    np.savetxt('../hsc_%d.txt' % nside, all)

    ##  Celestial bodies.   
    fig  = plt.figure(figsize=(10, 5))

    surveys       =  ['SO']
    masks, nside  =  prep_masks(surveys)
        
    def  update(month):
      pl.clf()

      ax   = fig.add_subplot(111, projection='mollweide', facecolor = 'black', alpha=0.1)

      plots(masks, surveys, clf=False, save=False)

      days = month * 30.

      RA, DEC    = get_plane()
      plot_mwd(ax, np.array(RA), np.array(DEC), org=0, color='r', alpha=1.0)

      RA, DEC    = get_plane(input=ephem.Ecliptic, output=ephem.Equatorial)
      plot_mwd(ax, np.array(RA), np.array(DEC), org=0, color='gold', alpha=1.0)

      ## ephem.Jupiter;  ['red', 'o']
      for body, props in zip([ephem.Sun, ephem.Moon], [['gold', '*'], ['gainsboro', 'o']]):
        body  = body(ephem.now()  + days)

        RA    = np.array([hours2degs(str(body.ra), dec=False)])
        DEC   = np.array([hours2degs(str(body.dec), dec=True)])

        plot_mwd(ax, RA, DEC, org=0, color=props[0], alpha=1.0, s=30, marker=props[1])

      pl.title(ephem.Date(ephem.now() + days))

    ## 12 * 1  
    ani        =  animation.FuncAnimation(fig, update, 3)
    writer     =  animation.writers['ffmpeg'](fps = 1)

    ani.save('demo.mp4', writer=writer, dpi=800)

    print('\n\nDone.\n\n')
