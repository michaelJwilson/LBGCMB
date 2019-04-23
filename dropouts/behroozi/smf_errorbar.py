import numpy             as np
import pylab             as pl
import matplotlib.pyplot as plt

from   schechterfn       import SchechterLfn
from   utils             import latexify


latexify(columns=1, equal=True, fontsize=12, ratio=None, ggplot=True, usetex=True)

##  Marchesini:  z=2.5                                                                                                                                                        
dat = np.loadtxt('dat/data-compilation/smf_ms/marchesini_z2.5.smf')

##  Columns:  Log10(stellar mass) (Msun), Log10(ND) (1 / Mpc^3 / dex), Err+ (dex), Err- (dex).                                                                                     
pl.errorbar(dat[:,0], dat[:,1], yerr=[dat[:,2], dat[:,3]], label = '2.5', markersize=2)

##  Marchesini:  z=3.5                                                                                                                                                        
dat = np.loadtxt('dat/data-compilation/smf_ms/marchesini_z3.5.smf')

##  Columns:  Log10(stellar mass) (Msun), Log10(ND) (1/Mpc^3/dex), Err+ (dex), Err- (dex)                                                                                     
pl.errorbar(dat[:,0], dat[:,1], yerr=[dat[:,2], dat[:,3]], label = '3.5', markersize=2)

##  Table 3 and Fig. (10) of https://arxiv.org/pdf/1507.05636.pdf 
##  Schechter fn. form of the Luminosity type.
'''
## z, log M* [Msun], alpha, phi* [Mpc^-3] after bottom rescaling. 
fits       = np.array([[4., 10.50, -1.55, 25.68], [5., 10.97, -1.70, 5.16], [6., 10.72, -1.91, 1.35]])
fits[:,1]  = 10. ** fits[:,1] 
fits[:,3] *= 1.e-5

for x in fits[:-1]:
  Ms       = np.logspace(7, 12, 50)
  Phis     = SchechterLfn(Ms, x[3], x[1], x[2])  ## [Mpc^-3]

  ##  Conversion from dM_* to dex, i.e. dM_* = dlog10 M_* x M_* x ln(10).
  ##  pl.loglog(Ms, Phis * np.log(10.) * Ms, label=x[0])
'''
##  Table 2 of https://arxiv.org/pdf/1507.05636.pdf. 
Ms = np.array([7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25])

## [Mpc^-3 / dex]
zs  = [4., 5., 6.]
Ps  = np.array([[-1.57, -1.77, -2.00, -2.22, -2.52, -2.91, -3.37,  -4.00,  -4.54],\
                [-1.47, -1.72, -2.01, -2.33, -2.68, -3.12, -3.47,  -4.12,  -4.88],\
                [-1.47, -1.81, -2.26, -2.65, -3.14, -3.69, -4.27, np.NaN, np.NaN]])

uEs = np.array([[0.21, 0.15, 0.13, 0.09, 0.09, 0.12, 0.09, 0.20, 0.34],\
                [0.24, 0.20, 0.16, 0.15, 0.07, 0.09, 0.16, 0.25, 0.40],\
                [0.35, 0.23, 0.21, 0.15, 0.12, 0.12, 0.38, np.NaN, np.NaN]])

dEs = np.array([[0.16, 0.14, 0.10, 0.09, 0.09, 0.05, 0.12, 0.25,   0.55],\
                [0.21, 0.20, 0.16, 0.10, 0.14, 0.11, 0.14, 0.38,   0.61],\
                [0.32, 0.28, 0.16, 0.15, 0.11, 0.13, 0.86, np.NaN, np.NaN]])

for kk, x in enumerate(Ps[:-1]):
  yerrs = [dEs[kk], uEs[kk]]

  pl.errorbar(Ms, x, yerr=yerrs, label=str(zs[kk]), markersize=2)

'''
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
'''

pl.ylim(-6.01, -0.99)

pl.xlabel(r'$\log_{10}(M_{*} \ / \ M_\odot)$')
pl.ylabel(r'$\log_{10}(\Phi \ / \ \rm{dex} \ / \ \rm{Mpc}^{-3})$')
 
##  handlelength=1, columnspacing=.2,
pl.legend(loc=1, ncol=2, frameon=False)

plt.tight_layout()

##  pl.show()
pl.savefig('plots/StellarMass_Schechter_errs.pdf')
