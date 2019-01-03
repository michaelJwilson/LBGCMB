import camb

import numpy as np
import pylab as pl

from   scipy              import optimize
from   scipy.interpolate  import interp1d
from   camb               import model, initialpower


params           = {'om_m':0.292, 'om_b':0.0468, 'om_L': 0.708, 'h_100':0.69, 'sig_8':0.819, 'Tcmb':2.786}  ## RunPB.
params['om_cm']  = params['om_m'] - params['om_b']

data             = np.loadtxt('EisensteinAndHu.dat')  ## See EisensteinAndHu.c for configuration choices.

ns               = 0.9688e+0
As               = 2.1870e-9

ks               = data[:,0] 
Tk               = data[:,1]
Tk_nu            = data[:,2]

PriPk            =   As * ks**ns ## Primordial, As = 1.0
EHpk             = PriPk * Tk**2.
EHpk_nu          = PriPk * Tk_nu**2.


## Load camb. 
pars  = camb.CAMBparams()

pars.set_accuracy(AccuracyBoost=2)

pars.InitPower.set_params(ns=0.9688, As=2.187e-9)

pars.set_cosmology(H0=100.*params['h_100'], ombh2 = params['om_b']*params['h_100']**2., omch2 = params['om_cm']*params['h_100']**2., mnu=0.00, omk=0.)

pars.set_matter_power(redshifts = [2.0], kmax=10.0) 

pars.NonLinear = model.NonLinear_none

linpk          = camb.get_results(pars)

sigma8         = linpk.get_sigma8()

k, z, linpk    = linpk.get_matter_power_spectrum(minkh=1e-4, npoints = 4000, maxkh=10.)
linpk          = linpk[0]


CAMBTK         = camb.get_transfer_functions(pars)

CAMBTK         = CAMBTK.get_matter_transfer_data()

CAMBTKk        = CAMBTK.transfer_data[0,:]
CAMBTK         = CAMBTK.transfer_data[6,:][:,0]
CAMBTK        /= CAMBTK[0]

'''
- Transfer_kh              = 1 (k/h)
- Transfer_cdm             = 2 (cdm)
- Transfer_b               = 3 (baryons)
- Transfer_g               = 4 (photons)
- Transfer_r               = 5 (massless neutrinos)
- Transfer_nu              = 6 (massive neutrinos)
- Transfer_tot             = 7 (total matter)
- Transfer_nonu            = 8 (total matter excluding neutrinos)
- Transfer_tot_de          = 9 (total including dark energy perturbations)
- Transfer_Weyl            = 10 (Weyl potential)
- Transfer_Newt_vel_cdm    = 11 (Newtonian CDM velocity)
- Transfer_Newt_vel_baryon = 12 (Newtonian baryon velocity)
- Transfer_vel_baryon_cdm  = 13 (relative baryon-cdm velocity)
'''

interp_linpk = interp1d(k, linpk,          bounds_error=False, fill_value=0.0)
interp_EHpk  = interp1d(ks, 2.9e14 * EHpk, bounds_error=False, fill_value=0.0)

kk           = np.arange(1e-4, 10.0, 1e-4)

def diff(A):
    interim  = (interp_linpk(kk[kk < 1e-2]) - A * interp_EHpk(kk[kk < 1e-2]))**2.
    
    result   = interim.sum()

    return result

result = optimize.minimize_scalar(diff, bounds=[0.9, 1.1])    

print result


interp_EHpk  = interp1d(ks, result.x * 2.9e14 * EHpk, bounds_error=False, fill_value=0.0)

pl.semilogx(kk, interp_linpk(kk)/interp_EHpk(kk) - 1, 'k-')

pl.xlim(1e-3, 1e0)

pl.title(r'$P(k) / P_{\rm{nw}}(k) - 1$')

pl.legend()

pl.savefig('wigglediff.pdf')

output = np.column_stack((kk, interp_linpk(kk) - interp_EHpk(kk)))
pl.savetxt("wigglediff.dat", output, header="k [h Mpc^-1]; P_L(k) - P_nw(k)")

output = np.column_stack((kk, interp_EHpk(kk)))
pl.savetxt("nowiggle.dat", output, header="k [h Mpc^-1]; P_nw(k)")

output = np.column_stack((kk, interp_linpk(kk)))
pl.savetxt("wiggle.dat", output, header="k [h Mpc^-1]; P_L(k)")
