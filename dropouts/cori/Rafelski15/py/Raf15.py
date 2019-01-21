import  os
import  numpy              as      np
import  pylab              as      pl

from    scipy.interpolate  import  interp1d
from    utils              import  latexify


def  Raf15_dNdz(droptype='u', field='UVUDF', depth='FULL', dz=0.1, no_lowz=False):
    root    = os.environ['LBGCMB'] + '/dropouts/cori/Rafelski15/dNdz/'

    ##  e.g. UVUDF_DEFAULTDEPTHS_udrops_dz_0.10.txt
    fname = '%s_%s_%sdrops_dz_%.2lf.txt' % (field, depth, droptype, dz)
    fpath = root + fname

    file  = np.loadtxt(fpath)
    midz  = file[:,0]
    dNdz  = file[:,1]

    if no_lowz:
        dNdz[midz < 1.0] = 0.0

    ##  np.sum(pz(zi)) * dz = 1.   
    ##  pz = interp1d(file[:,0], file[:,1], kind='linear', copy=True, fill_value=0.0, assume_sorted=False, bounds_error=False)
    
    return  midz, dNdz

def Raf15_pz(droptype='u', field='UVUDF', depth='FULL', dz=0.1, no_lowz=False):
    midz, dNdz = Raf15_dNdz(droptype=droptype, field=field, depth=depth, dz=dz, no_lowz=no_lowz)

    ngal       = np.sum(dNdz)
    ps         = dNdz / ngal / dz

    ##  print(field, np.sum(ps) * dz)
    return  midz, ps

def Raf15_intfrac(droptype='u', field='UVUDF', depth='FULL', dz=0.1, no_lowz=False):
    midz, dNdz = Raf15_dNdz(droptype=droptype, field=field, depth=depth, dz=dz, no_lowz=no_lowz)

    ngal       = np.sum(dNdz)
    nint       = np.sum(dNdz[midz < 1.6])

    return  ngal, nint


if __name__ == '__main__':
    import  pylab  as  pl


    print('\n\nWelcome to get dNdz.\n\n')

    dz     = 0.1
    all_ps = []

    latexify(columns=2, equal=False, fontsize=10, ratio=0.6, ggplot=True, usetex=True)

    ##  ['aegis', 'cosmos', 'goodsn', 'goodss', 'uds']     
    for depth, title in zip(['FULL', 'DEFAULTDEPTHS'], ['Full depth', 'Reduced depth']):
      pl.clf()

      for ttype, label in zip(['u', 'g', 'BzK', 'Euclid'], [r'$u$-dropouts', r'$g$-dropouts', r'$BzK$', 'All']):
        for field in ['UVUDF']:
          ngal, nint =  Raf15_intfrac(droptype=ttype, field=field, depth=depth, dz=dz) 
          midz, ps   =  Raf15_pz(droptype=ttype, field=field, depth=depth, dz=dz, no_lowz=True)

          if label is not 'All':
              label += ' (%d / %d)' % (nint, ngal)

          else:
              label += ' (%d)'      %  ngal

          pl.plot(midz, ps, '-', label=label, alpha=0.75)

          ## Collect for mean over fields:
          ## all_ps.append(ps)

          '''                                                                                                                                       
          all_ps  = np.c_[all_ps]                                                                                                                    
          mean_ps = np.mean(all_ps, axis=0)                                                                                                          
          std_ps  =  np.std(all_ps, axis=0)                                                                                                           
                                                                                                                                                
          pl.plot(midz, mean_ps, 'k-')                                                                                                                
          pl.errorbar(midz, mean_ps, std_ps, c='k')                                                                                                    
          '''

      pl.xlabel(r'$z$')
      pl.ylabel(r'$p(z)$')

      pl.title(title)

      pl.legend(ncol=2)

      pl.show()
      ## pl.savefig(os.environ['LBGCMB'] + '/dropouts/cori/Rafelski15/plots/pz_%s.pdf' % depth)

    print('\n\nDone.\n\n')
