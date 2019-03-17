@np.vectorize
def const_nz(z, ngal=1.e-4, zmin=3.0, zmax=4.0):
    if (zmin <= z) & (z <= zmax):
      return  ngal     ##  [(h^{-1} Mpc)^-3]                                                                                                                                      

    else:
      return  0.0
