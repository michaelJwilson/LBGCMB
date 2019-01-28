import  os
import  numpy              as      np
import  pylab              as      pl 

from    scipy.interpolate  import  interp1d


def get_contamination(m, zee=4, depth='W'):
    '''                                                                                                                                                     
    Get magnitude bins defined by GoldRush and corresponding contamination fraction.                                                                        
    '''

    root          = os.environ['LBGCMB']

    if zee == 3:
      path        = root + '/dropouts/dat/yoshida08.dat'

    elif zee in [4., 5.]:
      path        = root + '/dropouts/goldrush/cats/contamination/contamination_z%d%s.dat' % (zee, depth)

    else:
        raise  ValueError("\n\nContamination rate is not available for z = %.4lf.\n\n" % zee)

    data           = np.loadtxt(path)

    magbins        = data[:,0]
    contamination  = data[:,1]

    interp         = interp1d(magbins, contamination, kind='linear', axis=-1, fill_value=(contamination[0], 0.0), bounds_error = False)

    return  interp(m)

if __name__ == '__main__':
    ms = np.arange(20., 26.0, 0.1)
    
    for depth in ['W', 'D']:
      for zee in [4., 5.]:
        cs = get_contamination(ms, zee=zee, depth=depth)

        pl.plot(ms, cs, label='%s%.1lf' % (depth, zee))

    pl.xlabel(r'$z$')
    pl.ylabel(r'Contamination')

    pl.legend()
    pl.show()

    print('\n\nDone.\n\n')
