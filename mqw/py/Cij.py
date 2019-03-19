def Cij(Pk_interps, Llls, zmin, zmax):
    zee       = (zmin + zmax) / 2.

    chi       = comoving_distance(zee)
    dchi      = comoving_distance(zmax) - comoving_distance(zmin)

    ##  Extended Limber approximation.                                                                                                                         
    ks        = (Llls + 0.5) / chi                       ##  For the Phh evaluation in the integral, we take a zeff approx.                                    
                                                         ##  i.e. \int dz .... Phh(zeff).                                                                      

    ##  Effectively measures the non-linear matter power spectrum, but normalised down (by flux factor).                                                       
    return  Pmm(Pk_interps, ks, zee, 'nlinear') / chi ** 2. / dchi
