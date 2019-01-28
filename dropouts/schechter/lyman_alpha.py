import  numpy        as      np
import  pylab        as      pl

from    params       import  get_params
from    schechterfn  import  SchechterLfn


params    =   get_params()

##  Taken from Ly-a LF at z = 2.2 (1512.01854)
phi_star  =   6.32e-4                 ## [Mpc   ** -3.]
phi_star /=   params['h_100'] ** 3.   ## [Mpc/h ** -3.]  

print  phi_star * 1.e3

alpha     =  -1.75

LSolar    =  3.846e33                 ## [ergs / s]
Lstar     =  5.290e42                 ## [ergs / s]  

Lstar    /=  LSolar                   ## [L solar]

Ls        = 10. ** np.arange(6., 11., 0.01)           ## [L solar]
Phis      = SchechterLfn(Ls, phi_star, Lstar, alpha)  ## [h_100 / Mpc]^3

pl.loglog(Ls, Phis)

pl.xlabel(r'L [$L_\odot$]')

pl.savefig('plots/Lyman_alpha.pdf')
