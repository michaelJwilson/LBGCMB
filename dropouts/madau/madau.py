import  numpy  as      np


def lephare_madau(rwl, z):
    """
    Input:   Rest wavelength, redshift. 

    Output:  Return the madau extragalactic extinction (by neutral hydrogen) curve
             shipped with Le Phare.
    """

    from scipy  import interpolate

    if z > 8.0:
        raise ValueError("Le Phare Madau correction is not available for z > 8.0.")

    lo_extension  = np.arange(0.0, 18., 1.)
    lo_ones       = np.ones_like(lo_extension)

    hi_extension  = np.arange(1216., 250000., 1.)
    hi_ones       = np.ones_like(hi_extension)

    data          = np.loadtxt("LePhare/tau{0:02d}.out".format(np.int(10 * z)))

    data          = np.concatenate([np.column_stack((lo_extension, lo_ones * data[0,1])), data])
    data          = np.concatenate([data, np.column_stack((hi_extension, hi_ones))])
    
    interp        = interpolate.interp1d(data[:,0], data[:,1], bounds_error=False, fill_value=1.0)

    return interp(rwl)

def efftau_madau(rwl, z):
    """
    Madau 1995 extinction for a galaxy spectrum at redshift z, defined on a wavelenght grid wl [A].
    https://github.com/spacetelescope/synphot_refactor/issues/77

    This approximation is from Footnote 3 of Madau et al. 1995, ApJ, 441, 18. This is claimed accurate to 5%. 
    And the scatter in this factor (due to different lines of sight) is huge (attached see Madau's Fig. 3 top panel). 
    Also see Madau's Fig. 3 bottom panel for a redshifted version of the 'exact' prescription.
    """

    from  numpy  import  array, where, exp

    la  = 1215.67    ## Lyman alpha. Angstroms.                                                                                                         
    lb  = 1026.      ## Lyman beta.  Angstroms.                                                                                                        
    lg  = 973.
    ld  = 950.
    le  = 938.
    ll  = 912.       ## Lyman limit.

    ## Redshifted from restframe to (definitely) observed. 
    wl  = rwl*(1. + z)
    
    n   = len(wl)
    c   = array([3.6e-3, 1.7e-3, 1.2e-3, 9.3e-4])
    l   = array([    la,     lb,     lg,     ld])

    tau = np.zeros_like(wl)
    xe  = 1. + z

    ## Lyman series
    for i in range(len(l) ):
        indices       = where(wl <= l[i]*xe)           ## Note: no lower wavelength limit is correct as 'broadcasted' corrections. 
        tau[indices] += c[i]*(wl[indices]/l[i])**3.46

    ## Photoelectric absorption
    xc  = wl/ll
    xc3 = xc**3
    
    tau = where(wl <= ll*xe, tau + 0.25*xc3*(xe**.46 - xc**0.46) \
                                 + 9.4*xc**1.5*(xe**0.18 - xc**0.18) \
                                 - 0.7*xc3*(xc**(-1.32) - xe**(-1.32)) \
                                 - 0.023*(xe**1.68-xc**1.68), tau)
    
    '''
    min_tau     = tau.min()
    index       = np.where(tau == min_tau)
    tau[:index] = min_tau 
    
    tau         = where(tau < 0.0, 0.0, tau)
    '''

    return where(tau > 700., 0., exp(-tau))


if __name__ == "__main__":
    import  pylab as pl


    print("\n\nWelcome to Madau.\n\n")

    colors        = ['k', 'c', 'r', 'b', 'g', 'y', 'm']

    lo_extension  = np.arange(0.0, 18., 1.)
    lo_ones       = np.ones_like(lo_extension)

    hi_extension  = np.arange(1216., 10000., 1.)
    hi_ones       = np.ones_like(hi_extension)

    for c, z in zip(['k'], [6.5]):
    ## for c, z in zip(colors, [2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.]):
      ## Observed wavelength, i.e. redshifted - change source. 
      xs = np.arange(0., 15550., 1.)
      rs = xs/(1. + z)

      ## pl.plot(rs, efftau_madau(rs, z), label='z=%s' % z, color=c)
      ## pl.plot(xs, efftau_madau(rs, z), label='z=%s' % z)

      '''
      pl.axvline((1. + z)*la, c='k', ymin=0., ymax=1., lw=0.1)
      pl.axvline((1. + z)*lb, c='k', ymin=0., ymax=1., lw=0.1)
      pl.axvline((1. + z)*lg, c='k', ymin=0., ymax=1., lw=0.1)
      pl.axvline((1. + z)*ld, c='k', ymin=0., ymax=1., lw=0.1)
      pl.axvline((1. + z)*le, c='k', ymin=0., ymax=1., lw=0.1)
      pl.axvline((1. + z)*ll, c='k', ymin=0., ymax=1., lw=0.1)
      '''

      pl.plot(rs, lephare_madau(rs, z), color=c, linestyle='-', label='z=%s' % z)

    pl.xlim(0.0, 1300.)
    ## pl.xlim(0.0, 8550.)
    
    pl.ylim(-0.05, 1.05)

    ## pl.ylabel(r'$\tau_{\rm{eff}}$')
    pl.ylabel(r'$\rm{exp}(-\tau_{\rm{eff}})$')

    pl.xlabel(r'$\lambda_{\rm{rest}}$ [Angstroms]')
    ## pl.xlabel(r'$\lambda_{\rm{obs}}$ [Angstroms]')

    pl.legend(loc=2)

    pl.savefig('plots/madau.pdf')

    print("\n\nDone.\n\n")
    
