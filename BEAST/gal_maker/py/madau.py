import  os
import  numpy                  as     np
import  astropy.units          as     u

from    desispec.interpolation import resample_flux
from    magABsource            import magAB


la = 1215.67    ## Lyman alpha. Angstroms.                                                                                                                                                                                                                                          
lb = 1026.      ## Lyman beta.  Angstroms.                                                                                                                                                                                                                                          
lg = 973.
ld = 950.
le = 938.
ll = 912.       ## Lyman limit.   

def lephare_madau(rwl, z):
    '''
    Input:   Rest wavelength, redshift. 
    Output:  Return the madau extragalactic extinction (by neutral hydrogen) 
             curve shipped with Le Phare.
    '''

    if z > 8.0:
        raise ValueError('Le Phare Madau correction is not available for z > 8.0.')

    ##  http://www.cfht.hawaii.edu/~arnouts/LEPHARE/download.html
    root   = os.environ['BEAST']
    dat    = np.loadtxt(root + '/gal_maker/Madau/tau{0:02d}.out'.format(np.int(10 * z)))

    lext   = np.arange(0.0,    18., 1.)
    hext   = np.arange(1216., 1.e4, 1.)



    result = resample_flux(rwl, np.concatenate([lext, dat[:,0], hext]), np.concatenate([dat[0,1] * np.ones_like(lext), dat[:,1], np.ones_like(hext)]), extrapolate=True)

    return  result

def efftau_madau(rwl, z):
    '''
    Madau 1995 extinction for a galaxy spectrum at redshift z, defined on a wavelength grid wl [A].
    https://github.com/spacetelescope/synphot_refactor/issues/77

    Below is the Lyman-alpha prescription from the photo-z code BPZ.

    This approximation is from Footnote 3 of Madau et al. 1995, ApJ, 441, 18. This is claimed accurate to 5%. 
    And the scatter in this factor (due to different lines of sight) is huge (attached see Madau's Fig. 3 top panel). 
    Also see Madau's Fig. 3 bottom panel for a redshifted version of the 'exact' prescription.
    '''

    from  numpy  import  array, where, exp


    ##  Redshifted from restframe to (definitely) observed. 
    xe  = (1. + z)
    wl  = rwl * xe
    
    n   = len(wl)
    c   = array([3.6e-3, 1.7e-3, 1.2e-3, 9.3e-4])
    l   = array([    la,     lb,     lg,     ld])

    tau = np.zeros_like(wl)

    ##  Lyman series
    for i in range(len(l)):
        ##  Note: no lower wavelength limit is correct as 'broadcasted' corrections.
        indices       = where(wl <= l[i] * xe)            
        tau[indices] += c[i] * (wl[indices] / l[i])**3.46

    ##  Photoelectric absorption
    xc  = wl / ll
    xc3 = xc**3
    
    tau = where(wl <= ll*xe, tau + 0.25*xc3*(xe**.46 - xc**0.46) \
                                 + 9.4*xc**1.5*(xe**0.18 - xc**0.18) \
                                 - 0.7*xc3*(xc**(-1.32) - xe**(-1.32)) \
                                 - 0.023*(xe**1.68-xc**1.68), tau)
    
    return where(tau > 700., 0., exp(-tau))


if __name__ == "__main__":
    import  pylab as pl


    print("\n\nWelcome to Madau.\n\n")

    colors        = ['k', 'c', 'r', 'b', 'g', 'y', 'm', 'orange', 'm']

    waves         = np.arange(3.6e3, 1.e4, 0.1)
    
    uwaves        = waves * u.AA
    vs            = uwaves.to(u.Hz, equivalencies=u.spectral()).value 
    Fv            = 1. * np.ones_like(vs) ##  [ergs/s/cm2/Hz]

    for c, z in zip(colors, np.arange(0.0, 9.0, 1.0)):
      rs = waves / (1. + z)

      pl.axvline((1. + z) * la, c='k', ymin=0., ymax=1., lw=0.1, alpha=0.2)
      pl.axvline((1. + z) * lb, c='k', ymin=0., ymax=1., lw=0.1, alpha=0.2)
      pl.axvline((1. + z) * lg, c='k', ymin=0., ymax=1., lw=0.1, alpha=0.2)
      pl.axvline((1. + z) * ld, c='k', ymin=0., ymax=1., lw=0.1, alpha=0.2)
      pl.axvline((1. + z) * le, c='k', ymin=0., ymax=1., lw=0.1, alpha=0.2)
      pl.axvline((1. + z) * ll, c='k', ymin=0., ymax=1., lw=0.1, alpha=0.2)

      pl.plot(waves, lephare_madau(rs, z) * Fv, c, label='z=%s' % z)
      pl.plot(waves,  efftau_madau(rs, z) * Fv, 'k--')

    pl.xlim(3.6e3, 1.e4)
    pl.ylim(1.e-7, 2.00)

    pl.xlabel(r'$\lambda$ [Angstroms]')
    pl.ylabel(r'$F_\nu$')

    pl.yscale('log')

    pl.legend(loc=4, ncol=3)
    pl.show()

    print("\n\nDone.\n\n")
    
