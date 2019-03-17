def check_nP(Pk_interps, fsky=0.5, nz=const_nz):
    zs =   np.arange(0.0, 6.0, 0.5)
    ks = np.logspace( -2,   0,  20)

    for z in zs:
      pl.semilogx(ks, nP(z, 0.0, ks, Pk_interps, fsky=fsky, nz=nz), label=str(z))

    pl.xlabel(r'$k$')
    pl.ylabel(r'$nP/(1 + nP)$')

    pl.legend()
    pl.show()

def check_dropnz(nz):
    zs = np.arange(0., 10., 0.01)

    pl.plot(zs, drop_nz(zs), 'k-')

    pl.xlabel(r'$z$')
    pl.ylabel(r'$\bar n(z)$')

    plt.tight_layout()
    pl.show()

def check_vol(fsky=0.5, fkp_weighted=False, nbar=1.e-3, P0=5.e3):
    ##  Simple volume integral as a sanity check.                                                                                                                                
    zmin    = 0.6
    zmax    = 1.2

    zranges = [[zmin, zmax]]
    args    = (fsky, fkp_weighted, nbar, P0)

    ##  vol_integrand(z, fsky=0.5, fkp_weighted=False, nbar=1.e-3, P0=5.e3)                                                                                                   
    result  = nquad(vol_integrand, zranges, args=args, full_output=True)

    args    = (fsky, True, nbar, P0)
    fkp_wt  = nquad(vol_integrand, zranges, args=args, full_output=True)

    print('Vol: %.4lf [(h^-1 Gpc)^3], FKP Vol:  %.4lf [(h^-1 Gpc)^3], compared to Astropy: %.4lf [(h^-1 Gpc)^3]' % (result[0] / 1.e9, fkp_wt[0] / 1.e9,\
                                                                                                                    fsky * (cosmo.comoving_volume(zmax).value - cosmo.comoving_volume(zmin).value) * cparams['h_100'] ** 3. / 1.e9))

def _vcheck_vol(fsky=0.5, fkp_weighted=True, nbar=1.e-3, P0=5.e3):
    zmin    =  0.6
    zmax    =  1.2

    zranges =  [[zmin, zmax]]

    integ   =  vegas.Integrator(zranges)
    args    = (fsky, fkp_weighted, nbar, P0)

    result  =  integ(lambda x: _vvol_integrand(x, args), nitn=10, neval=10000)

    print(result)
    print('\n\nAstropy: %.4lf [(h^-1 Gpc)^3]' % (fsky * (cosmo.comoving_volume(zmax).value - cosmo.comoving_volume(zmin).value) * cparams['h_100'] ** 3.))
