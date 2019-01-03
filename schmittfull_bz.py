import  os
import  numpy              as      np
import  pylab              as      pl
import  pandas             as      pd

from    growth_rate        import  growth_factor
from    collections        import  OrderedDict
from    scipy.interpolate  import  interp1d


def sdss(zee):
    zee     = np.asarray(zee)
    result  = 1.0 + (zee - 0.1)

    ## result[zee < 0.1] = 1.0

    return  result

def boss_lrg(zee, force=False):
    zee = np.asarray(zee)
    aa  = 1. / (1. + zee)

    if force:
        return  1.70 / growth_factor(aa)

    if np.any(zee > 1.0):
      raise ValueError('\n\nCalling b(z) outside range of BOSS LRGs (z < 1.0).\n\n')

    return  1.70 / growth_factor(aa)

def desi_bgs(zee, force=False):
    zee = np.asarray(zee)
    aa  = 1. / (1. + zee)

    if force:
        return  1.34 / growth_factor(aa)

    if np.any(zee > 0.5):
      raise ValueError('\n\nCalling b(z) outside range of BGS (z < 0.5).\n\n')

    return  1.34 / growth_factor(aa)

def desi_lrg(zee, force=False):
    zee = np.asarray(zee)
    aa  = 1. / (1. + zee)

    if force:
        return  1.70 / growth_factor(aa)

    if np.any(zee > 1.3):
      raise ValueError('\n\nCalling b(z) outside range of LRGs (z < 1.3).\n\n')

    return  1.70 / growth_factor(aa)

def desi_elg(zee, force=False):
    zee = np.asarray(zee)
    aa  = 1. / (1. + zee)
    
    if force:
        return  0.84 / growth_factor(aa)

    if np.any(zee > 1.8):
      raise ValueError('\n\nCalling b(z) outside range of ELGs (z < 1.8).\n\n')

    return  0.84 / growth_factor(aa)

def desi_qso(zee, force=False):
    zee = np.asarray(zee)
    aa  = 1. / (1. + zee)

    if force:
        return  1.20 / growth_factor(aa)

    if np.any(zee > 5.0):
      raise  ValueError('\n\nCalling b(z) outside range of QSOs (z < 5.0).\n\n')
    
    return  1.20 / growth_factor(aa)

def desi_eff(zee):
    root = os.environ['LBGCMB']

    data = np.loadtxt(root + '/cadence/dat/effective_desi.dat')
    ibz  = interp1d(data[:,0], data[:,2], kind='linear', bounds_error=False, fill_value=0.0)

    zee  = np.asarray(zee)

    return  ibz(zee)

def get_allbz():
    '''
    0.0   0.5        22917     SDSS1
    0.5   0.8        22917     SDSS2
    0.0   0.5        685       DESIBGS
    0.0   0.9        139       BOSSLRGs
    0.6   1.2        279       DESILRGs
    0.6   0.8        250       DESIELGs1
    0.8   1.2        929       DESIELGs2
    0.6   1.9        100       DESIQSOs
    '''

    result               = OrderedDict()

    ##  FIX:  Cheap Hack. 
    result['SDSS1']      = lambda z: sdss(z, force=True)
    result['SDSS2']      = lambda z: sdss(z, force=True)
    result['DESIBGS']    = lambda z: desi_bgs(z, force=True)
    result['BOSSLRGs']   = lambda z: boss_lrg(z, force=True)
    result['DESILRGs']   = lambda z: desi_lrg(z, force=True)
    result['DESIELGs1']  = lambda z: desi_elg(z, force=True)
    result['DESIELGs2']  = lambda z: desi_elg(z, force=True)
    result['DESIQSOs']   = lambda z: desi_qso(z, force=True)
    
    return  result


if __name__ == '__main__':
    print('\n\nWelcome to Schmittfull & Seljak linear b(z).')

    zs = np.arange(0.2, 5.0, 0.01)

    pl.plot(zs,     sdss(zs), label='SDSS')
    pl.plot(zs, boss_lrg(zs), label='BOSS LRG')
    pl.plot(zs, desi_bgs(zs), label='DESI BGS')
    pl.plot(zs, desi_lrg(zs), label='DESI LRG')
    pl.plot(zs, desi_elg(zs), label='DESI ELG')
    pl.plot(zs, desi_qso(zs), label='DESI QSO')
    
    pl.legend(ncol=2)

    pl.xlabel(r'$z$')
    pl.ylabel(r'$b(z)$')
    
    pl.savefig('plots/ss17_biasz.pdf')
    
    print('\n\nDone.\n\n')
