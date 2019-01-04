import  os
import  numpy              as      np
import  pylab              as      pl

from    scipy.interpolate  import  interp1d
from    utils              import  latexify

def  get_Rafelski15_dNdz(droptype='u', field='aegis', depth='FULL', dz=0.1):
    root    = os.environ['LBGCMB'] + '/dropouts/cori/Rafelski15/dNdz/'

    ## e.g. cosmos_FULL_udrops_dz_0.10.txt
    fname = '%s_%s_%sdrops_dz_%.2lf.txt' % (field, depth, droptype, dz)
    fpath = root + fname

    file  = np.loadtxt(fpath)
    midz  = file[:,0]
    dNdz  = file[:,1]

    ##  np.sum(pz(zi)) * dz = 1.   
    ## pz    = interp1d(file[:,0], file[:,1], kind='linear', copy=True, fill_value=0.0, assume_sorted=False, bounds_error=False)
    
    return  midz, dNdz

def  get_Rafelski15_pz(droptype='u', field='aegis', depth='FULL', dz=0.1):
    midz, dNdz = get_Rafelski15_dNdz(droptype=droptype, field=field, depth=depth, dz=dz)

    ngal       = np.sum(dNdz)
    ps         = dNdz / ngal / dz

    ## print(field, np.sum(ps) * dz)
    return midz, ps


if __name__ == '__main__':
    import  pylab  as  pl


    print('\n\nWelcome to get dNdz.\n\n')

    dz     = 0.1
    fields = ['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']

    all_ps = []

    latexify(columns=2, equal=False, fontsize=10, ratio=None, ggplot=True, usetex=True)

    for field in fields:
      midz, ps = get_Rafelski15_pz(droptype='u', field=field, depth='FULL', dz=dz)
      all_ps.append(ps)

      pl.plot(midz, ps, '-', label=field.upper(), alpha=0.5)

    all_ps  = np.c_[all_ps]  
    mean_ps = np.mean(all_ps, axis=0)
    std_ps  =  np.std(all_ps, axis=0)

    pl.plot(midz, mean_ps, 'k-')
    pl.errorbar(midz, mean_ps, std_ps, c='k')

    pl.legend(ncol=2)
    pl.show()

    print('\n\nDone.\n\n')
