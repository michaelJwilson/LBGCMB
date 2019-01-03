import  numpy              as      np

from    specs              import  samplestats         as  gsample_stats
from    ilim               import  get_nbar_nocontam
from    bz                 import  get_dropoutbz
from    completeness       import  get_dropoutpz       as  get_gdropoutpz
from    scipy.interpolate  import  interp1d
from    astropy.table      import  Table, Column


band         =  'g'
depth        =  25.5
stats        =  gsample_stats()

## Assume no contamination in the UD field.
## stats     =  get_nbar_nocontam(band, depth='UD', printit=False)
## nbar      =  stats[band]['nbar_nointerlopers']

nbar         =  stats[band]['nbar']
peakz        =  stats[band]['z']

bz           =  get_dropoutbz(m=depth)

zee, pzee    =  get_gdropoutpz()
pz           =  interp1d(zee, pzee, kind='linear', bounds_error=False, fill_value=0.0, assume_sorted=False)


if __name__ == '__main__':
    print('\n\nWelcome Simo!\n\n')

    zs            =  np.linspace(0.1, 5.0, 20)
    ps            =  pz(zs)
    ps           /=  np.sum(ps)

    ns            =  nbar * ps
    bs            =  bz(zs)

    ind           =  np.where(ps == 0.0)
    bs[ind]       =  0.0

    np.savetxt('./dat/BEAST_dropband_%s_depth_%.1lf.txt' % (band, depth), np.c_[zs, ps, ns, bs, nbar * np.ones_like(zs)], fmt='%.6lf \t', header='zs\t\t p(z)\t\t n(z) [1/deg2]\t b(z) \t\t n(z < oo)')

    print(peakz, nbar, np.sum(ns))

    print('\n\nDone.\n\n')
    
