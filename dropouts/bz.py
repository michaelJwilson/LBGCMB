import  os
import  numpy              as      np

from    scipy.interpolate  import  interp1d
from    growth_rate        import  growth_factor


def get_dropoutbz(m=24.5):
    ##  z, b(z).
    data = np.loadtxt(os.environ['LBGCMB'] + '/dropouts/dat/dropout_bz.dat')

    if m == 24.5:
      return  interp1d(data[:-1,0], data[:-1,1], kind='quadratic', bounds_error=True, fill_value=0.0)

    elif m == 25.0:
      return  interp1d(data[:,0], data[:,2], kind='cubic', bounds_error=True, fill_value=0.0)

    elif m == 25.5:
      return  interp1d(data[:-1,0], data[:-1,3], kind='quadratic', bounds_error=True, fill_value=0.0)  

    else:
        raise  ValueError("Requested b(z) is not available.")
  

if __name__ == '__main__':
    import  pylab  as  pl


    print('\n\nWelcome to dropout b(z).\n\n')

    ##  Schmittfull and Seljak form.                                                                                                                  
    zs = np.arange(2.0, 6.0, 0.01)
    pl.plot(zs, (1. + zs), 'k-', label='$(1+z)$', alpha=0.6)

    ##  m = 24.5
    bz = get_dropoutbz()

    zs = np.arange(2.0, 5.0, 0.01)
    bs = bz(zs)

    pl.plot(zs, bs, label='m=24.5')

    ##  m = 25.0
    bz = get_dropoutbz(m=25.0)
    zs = np.arange(2.0, 6.0, 0.01)
    bs = bz(zs)

    pl.plot(zs, bs, label='m=25.0')

    ##  m = 25.5
    bz = get_dropoutbz(m=25.5)
    zs = np.arange(2.0, 5.0, 0.01)
    bs = bz(zs)

    pl.plot(zs, bs, label='m=25.5')

    pl.xlabel(r'$z$')
    pl.ylabel(r'$b(z)$')
    
    pl.legend()
    
    pl.show()
    ##  pl.savefig('plots/dropout_bz.pdf')

    print('\n\nDone.\n\n')
