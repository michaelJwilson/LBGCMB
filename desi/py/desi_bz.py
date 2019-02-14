import os 
import numpy as np

from   growth_rate import growth_factor


def lrg_bz(zz):
    aa = 1. / (1. + zz)
    Dp = growth_factor(aa)

    return 1.7 / Dp

def elg_bz(zz):
    aa = 1. / (1. + zz)
    Dp = growth_factor(aa)

    return 0.84 / Dp

def bgs_bz(zz):
    aa = 1. / (1. + zz)
    Dp = growth_factor(aa)

    return 1.0 / Dp

def qso_bz(zz):
    aa = 1. / (1. + zz)
    Dp = growth_factor(aa)

    return 1.2 / Dp


if __name__ == '__main__':
    import pylab as pl


    print('\n\nWelcome to DESI b(z). \n\n')

    zs = np.arange(0.3, 0.6, 0.1)
    bs = bgs_bz(zs)

    pl.plot(zs, bs)
    pl.show()
    
    print('\n\nDone.\n\n')
