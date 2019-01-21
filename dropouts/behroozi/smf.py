import numpy             as np
import pylab             as pl
import matplotlib.pyplot as plt

from   schechterfn       import SchechterLfn
from   utils             import latexify


latexify(columns=1, equal=True, fontsize=10, ratio=None, ggplot=True, usetex=True)

##  Table 3 and Fig. (10) of https://arxiv.org/pdf/1507.05636.pdf 
##  Schechter fn. form of the Luminosity type.

## z, log M* [Msun], alpha, phi* [Mpc^-3] after bottom rescaling. 
fits       = np.array([[4., 10.50, -1.55, 25.68], [5., 10.97, -1.70, 5.16], [6., 10.72, -1.91, 1.35]])
fits[:,1]  = 10. ** fits[:,1] 
fits[:,3] *= 1.e-5

for x in fits:
  Ms       = np.logspace(7, 12, 50)
  Phis     = SchechterLfn(Ms, x[3], x[1], x[2])  ## [Mpc^-3]

  ##  Conversion from dM_* to dex, i.e. dM_* = dlog10 M_* x M_* x ln(10).
  pl.loglog(Ms, Phis * np.log(10.) * Ms, label=x[0])

##  Table 2 of https://arxiv.org/pdf/1507.05636.pdf. 
Ms = np.array([7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25])

## [Mpc^-3 / dex]
zs = [4., 5., 6.]
Ps = np.array([[-1.57, -1.77, -2.00, -2.22, -2.52, -2.91, -3.37,  -4.00,  -4.54],\
               [-1.47, -1.72, -2.01, -2.33, -2.68, -3.12, -3.47,  -4.12,  -4.88],\
               [-1.47, -1.81, -2.26, -2.65, -3.14, -3.69, -4.27, np.NaN, np.NaN]])

for kk, x in enumerate(Ps):
  pl.loglog(10. ** Ms, 10. ** x, '^', label=zs[kk], markersize=2)

##  Table 1 of https://arxiv.org/pdf/1303.4409.pdf
##  Schechter fn. form of the Luminosity type.
##  [z, log(M^*_star), Phi^* [1e-4 Mpc^{-3}], alpha, logMlim_star]
Muzzin       = np.array([[2.75, 11.03, 1.93, -1.01, 10.76], [3.5, 11.49, 0.09, -1.45, 10.94]])
Muzzin[:,2] *= 1.e-4  ## [Mpc^{-3}]

def MuzzinSchechter(StellarMass, PhiStar, Mstar, alpha):
  ##  Stellar Mass [Msun]; Mstar = np.log10(M^*_star / Msolar).
  M = np.log10(StellarMass)

  ##  eqn. (4) of https://arxiv.org/pdf/1303.4409.pdf
  return  np.log(10.) * PhiStar * 10. ** ((1. + alpha) * (M - Mstar)) * np.exp(-10. * (M - Mstar))

StellarMass = np.logspace(Muzzin[1,4], 12, 50)
pl.loglog(StellarMass, MuzzinSchechter(StellarMass, Muzzin[1,2], Muzzin[1,1], Muzzin[1,3]), 'c-', label='Muzzin')

pl.ylim(1.e-5, 1.e-1)

pl.xlabel(r'$M_{*} \ [M_\odot]$')
pl.ylabel(r'$\Phi \ [(\rm{Mpc})^{-3} \rm{dex}^{-1}]$')
 
pl.legend(loc=3, ncol=2)
plt.tight_layout()
pl.savefig('plots/StellarMass_Schechter.pdf')

pl.clf()

##  Cumulative \bar n.
dx = 0.1

xs = np.arange(5, 15., dx)
Ms = 10. ** xs

log10Mmins = np.arange(9., 11., 0.1)

for x in fits:
  Phis     = SchechterLfn(Ms, x[3], x[1], x[2])  ## [Mpc^-3]                                                                                               
  result   = []

  for log10Mmin in log10Mmins:
    Mmin   = 10.** log10Mmin
    nbar   = np.sum(Phis[Ms >= Mmin] * Ms[Ms >= Mmin] * np.log(10.) * dx)

    result.append(nbar)

  result = np.array(result)

  ##  Conversion from dM_* to dex, i.e. dM_* = dlog10 M_* x M_* x ln(10).                                                                                
  pl.loglog(10. ** log10Mmins, result, label=x[0])

pl.xlabel(r'$M_* \ [M_\odot]$')
pl.ylabel(r'$\bar n \ [(\rm{Mpc})^{-3}]$')

pl.legend()
plt.tight_layout()
pl.savefig('plots/StellarMass_nbar.pdf')
