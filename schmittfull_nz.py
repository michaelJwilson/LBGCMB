import  pylab              as      pl
import  numpy              as      np
import  pandas             as      pd

from    scipy.interpolate  import  interp1d
from    scipy.integrate    import  simps
from    collections        import  OrderedDict 
from    schmittfull_bz     import  sdss, boss_lrg, desi_bgs, desi_lrg, desi_elg, desi_qso, get_allbz
from    whitebook_pz       import  const_pz


def ss_pz():
    ##  Schmittfull and Seljak (2017).
    data      = np.loadtxt('dat/schmittfull_nbar.dat')

    zs        = data[:,0]
    dndz      = data[:,1]
    
    nz_interp = interp1d(zs, dndz, kind='nearest', bounds_error=False, fill_value=0.0, assume_sorted=False) 
    
    zs        = np.arange(0., 9., 1.e-6)
    dndz      = nz_interp(zs)

    norm      = simps(dndz, dx = zs[1] - zs[0], axis=0)  ## [deg^-2]

    pz        = interp1d(zs, dndz / norm, kind='nearest', bounds_error=False, fill_value=0.0, assume_sorted=False)

    return  pz, 3600. * norm

def load_samples():
    return  pd.read_csv("dat/ss17_tracers.dat", names = ['zmin', 'zmax', 'N', 'Surveys'], sep = '\s+', comment='#')

def get_ss17_samples(nolsst=True):
    dframe    =  load_samples()
    surveys   =  list(dframe['Surveys'].values)

    zlos      =  dframe['zmin'].values.astype('float')
    zhis      =  dframe['zmax'].values.astype('float')

    all_bs    =  get_allbz()

    ns        =  dframe['N'].values.astype('float')
    bs        =  [all_bs[x] for x in surveys]

    ps        =  []

    for x in zip(zlos, zhis):
      pz      =   lambda z:  const_pz(z, zlo=x[0], zhi=x[1])  
      ps     +=  [np.vectorize(pz)]

    if not nolsst:
      ##  Add LSST sample.   
      pz, norm  = ss_pz()

      ns       += [norm]
      ps       += [pz]
      bs       += [lambda z:  (1. + z)]
      surveys  += ['LSST']

    print('\n\nAvailable surveys: ' + ''.join('  %s;' % s for s in surveys) + '\n\n')

    return  ns, ps, bs, surveys


if __name__ == '__main__':
    print('\n\nWelcome to Schmittfull <n>.\n\n')

    ## -- LSST -- ##
    tomos        = OrderedDict()

    tomos['1']   = [0.0, 0.5, 9.30e8]
    tomos['2']   = [0.5, 1.0, 1.55e9]
    tomos['3']   = [1.0, 2.0, 1.40e9]
    tomos['4']   = [2.0, 3.0, 2.40e8]
    tomos['5']   = [3.0, 4.0, 9.40e7]
    tomos['6']   = [4.0, 7.0, 4.30e7]

    pz, norm     = ss_pz()

    zs           = np.arange(0., 9., 1.e-6)
    ps           = pz(zs)

    survey_area  = 1.8e4
    
    for ii, key in enumerate(tomos):
      zlo        = tomos[key][0]
      zhi        = tomos[key][1]

      condition  = (zlo <= zs) & (zs <= zhi)

      result     = norm * simps(ps[condition], dx = zs[1] - zs[0], axis=0)

      print('%.1lf < z < %.1lf:  % 4.1lfM \t\t (% 4.1lfM)' % (zlo, zhi, result * survey_area / 1.e6, tomos[key][2] / 1.e6))
    
    print('\nTotal galaxies:  %.4lfM (%.4lfk)' % (survey_area * norm / 1.e6, norm / 1.e3))

    ## -- and the rest -- ##                                                                                       
    ## dframe       =  load_samples()
    ## surveys      =  list(dframe['Surveys'].values)
 
    ns, ps, bs, ss  =  get_ss17_samples(nolsst=True)

    '''
    samples         =  zip(bs, ps, ns)
    ## samples      =  samples[3:4]

    for p in ps:
      print p(0.5) 

    for [b, p, n] in samples:
      print p(0.5)
    '''

    print('\n\nDone.\n\n')
