import  os
import  numpy              as      np
import  pylab              as      pl

from    scipy.interpolate  import  interp1d


def  get_uvudf_dndz(band='g', field='UVUDF', mag=''):
    root    = os.environ['LBGCMB'] + '/dropouts/edison/dat/dNdz/'

    if mag is not '':
        mag = '_' + mag

    fname = band + field + mag + '.txt'
    fpath = root + fname

    file  = np.loadtxt(fpath)

    dz    = file[1,0] - file[0,0]

    ##  np.sum(pz(zi)) * dz = 1.   
    pz    = interp1d(file[:,0], file[:,1], kind='linear', copy=True, fill_value=0.0, assume_sorted=False, bounds_error=False)
    
    return  pz


if __name__ == '__main__':
    import  pylab  as  pl


    print('\n\nWelcome to get dNdz.\n\n')

    zs = np.arange(0.0, 10.0, 0.1)
    pz = get_uvudf_dndz(band='g', field='UVUDF', mag='')
    
    ps = pz(zs)

    pl.plot(zs, ps, 'k-')

    pl.show()

    print('\n\nDone.\n\n')
