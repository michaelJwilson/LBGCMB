#!/usr/bin/env python3
#
#  Fits to the stellar-mass-halo-mass relation from Ishikawa+17
#  https://xxx.lanl.gov/pdf/1612.06869v2
#  Eqs. (20, 21) and Table 3.
#  This uses the fitting function form of Behroozi+18, which I
#  reverse engineered from the "gen_smhm.py" script in the EDR.

import  pylab              as      pl
import  numpy              as      np
import  matplotlib         as      mpl
import  matplotlib.pyplot  as      plt

from    collections        import  OrderedDict
from    utils              import  latexify
from    params             import  get_params


latexify(fig_width=None, fig_height=None, columns=1, equal=True)

params = get_params()

def smhm_i17(Mh, drop_type):
    '''
    SMHM relation from Ishikawa+17 (arxiv:1612.06869; Eqs. 20,21 & Tab.3).
    Returns Mstar / Mh given Mh in Msun (not Msun/h).
    '''

    if drop_type == 'u':
        aa       = 0.25
        M1,eps   = 10**12.10, 2.222e-2

    elif drop_type == 'g':
        aa         = 0.20
        M1, eps    = 10**11.99, 2.248e-2

    elif drop_type == 'r':
        aa         = 0.17
        M1,eps     = 10**11.77, 2.091e-2

    else:
        raise ValueError('Unknown drop_type:  %s' % str(drop_type))

    ##  The defaults from Behroozi+18.
    lgMh  = np.log10(Mh)
    lna   = np.log(aa)
    lgM1  = 11.96766 + 2.230701*(aa-1)-2.358881*lna-0.3736004*(1./aa-1)
    rat1  = -1.357278 + 0.1385258*(aa-1)+0.2295319*lna+0.1565769*(1./aa-1)+lgM1
    alpha = 2.025043-1.364907*(aa-1)+1.174364*lna+0.1665663*(1./aa-1)
    beta  = 0.5204741-0.1346026*(aa-1)-0.1605413*(1./aa-1)
    gamma = 10**(-0.7288397-1.763953*(aa-1)-0.6393729*(1./aa-1))
    delta = 0.3512605

    ##  Overwrite two parameters with the results of Ishikawa+17, expressing
    ##  the normalization in terms of "eps".
    lgM1  = np.log10(M1)
    rat1  = np.log10(eps * M1) + np.log10(2.0) - gamma

    ##  and evaluate lg(Mstar/Mhalo):
    dM    = lgMh - lgM1
    f     = lambda dm: rat1 - np.log10( 10.**(-alpha * dm) + 10.**(-beta * dm) ) + gamma * np.exp(-0.5 * dm**2 / delta**2)

    return  f(dM) - lgMh


if __name__ == '__main__':
    print('\n\nWelcome to the stellar mass calculator.\n\n')
    
    Masses  = np.arange(11., 13.5, 0.01)
    results = OrderedDict()
    '''
    for drop_type, label in zip(['u', 'g', 'r'], [r'$u$', r'$g$', r'$r$']):
        results[drop_type] = []
        
        for lgMh in Masses:
            result = 10. ** smhm_i17(10. ** lgMh, drop_type)

            results[drop_type].append(result) 

            print(drop_type + ' {:.2f} {:8.5f}'.format(lgMh, result))

        Masses             = 10. ** Masses
        results[drop_type] = np.array(results[drop_type])

        pl.loglog(Masses, results[drop_type], label=label + '-dropouts')
    ''' 
   
    ##  SM / HM relation for dropouts from Ishikawa+17.
    ##  Quoting log10(Mstar / Mhalo) with Mhalo in Msun / h.
    data = np.loadtxt('../dat/stellar/smhm.txt')

    mh   = data[:,0]

    ud   = data[:,1]
    gd   = data[:,2]
    rd   = data[:,3]

    print mh
    print ud
    print gd
    print rd
    
    pl.loglog(mh, (mh / params['h_100']) * 10. ** ud, label=r'$u$-dropouts')
    pl.loglog(mh, (mh / params['h_100']) * 10. ** gd, label=r'$g$-dropouts')
    pl.loglog(mh, (mh / params['h_100']) * 10. ** rd, label=r'$r$-dropouts')
    
    pl.xlabel(r'$M_h \ [M_\odot / h]$')
    pl.ylabel(r'$M_\star \ [M_\odot]$')
    
    ## pl.ylim(5.e-4, 5.e-2)

    pl.legend()
    pl.savefig('../plots/stellar_mass.pdf', bbox_inches='tight')

    print('\n\nDone.\n\n')
