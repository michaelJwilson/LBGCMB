import numpy as np


def vipers_nz(z, A=3.103, z0=0.191, alpha=8.603, beta=1.448, ngal=5e3):
    ##  ngal is the expected number of galaxies per sq. deg.                                                                                                          
    ##  Eqn. (2) of https://arxiv.org/pdf/1303.2622.pdf                                                                                                               
    result  =  A * ((z / z0) ** alpha) * np.exp(-(z / z0) ** beta)

    ##  Normalisation.                                                                                                                                                
    dz      =  0.01
    zz      =  np.arange(0.0, 2.0, dz)
    norm    =  dz * np.sum(A * ((zz / z0) ** alpha) * np.exp(-(zz / z0) ** beta))

    return  result / norm

def plot_vipers(ngal=5.e3):
    ##  Test VIPERS N(z)                                                                                                                                                         
    zs   = np.arange(0.0, 1.3, 0.01)
    Ns   = vipers_nz(zs, A=3.103, z0=0.191, alpha=8.603, beta=1.448, ngal=5.e3)

    pl.plot(zs, Ns, label=ngal)
    pl.xlabel(r'$z$')
    pl.ylabel(r'$dN/dz$')
    pl.legend()
    plt.tight_layout()
    pl.show()
