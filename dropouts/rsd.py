import  numpy      as  np
import  pylab      as  pl

from    bz         import  get_dropoutbz
from    prep_camb  import  CAMB
from    pmh        import  Pmm, get_PkInterps, linz_bz
from    params     import  get_params


if __name__ == '__main__':
    z           = 4.0
    bz          = get_dropoutbz(m=25.0)

    params      = get_params()

    nbar        = np.array([6.e-4])
    nbar       /= params['h_100'] ** 3.

    cambx       =  CAMB()
    Pk_interps  =  get_PkInterps(cambx)

    ks          = np.logspace(-2., 0., 100.)
    Ps          = Pmm(Pk_interps, ks, z) * bz(z) * bz(z)

    pl.loglog(ks, Ps, label='z=3')
    pl.axhline(y=1./nbar, xmin=0., xmax=1.)
    
    pl.legend()

    pl.show()

    print('\n\nDone.\n\n')
