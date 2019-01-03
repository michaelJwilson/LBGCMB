#!/usr/bin/env python
#
# Python script to evaluate the number density as a function of
# apparent magnitude for U-band dropout galaxies at z~3.
#
from __future__ import print_function, division

import numpy as np


def print_nbar():
    """   
    Prints a table of lg10(nbar), in Mpc^{-3}, vs. apparent magnitude
    in the \mathcal{R} band for z~3 (2.7<z<3.4) galaxies.
    Reddy et al. 2008, ApJS, 175, 48, Table 7
    Phi* is in units of [h70/Mpc]^3/mag.
    K-correction from Eq. (14) of Reddy++ and preceding text.
    """
    print("# Number density of u-band dropout galaxies at z~3.")
    print("# Density vs mathcal{R}-band apparent magnitude.")
    print("# Assuming RunPB cosmology, units Mpc/h.")

    aa   = 0.25       # z=3

    chi  = 4483.6222  # Mpc/h for OmM=0.292, z=3.

    # h_100 in the RunPB cosmology.
    hub  = 0.69

    # the RunPB cosmology.
    dmod = 25. + 5.*np.log10(chi/aa/hub)

    print("# Distance modulus to z=%.1f is %6.2f."%(1./aa-1.,dmod))

    kcorr = 2.5*np.log10(1.0/aa)

    phis, mstar, alpha = 1.66e-3, -20.84, -1.57

    print("# %10s %10s %10s %10s" % ("mathcal{R}","M_AB(1700)","L/L*","lg[nbar]"))

    for maglim in [23.5, 24.0, 24.5, 25.0]:
        Mabs = maglim - dmod + kcorr

        L    = 10.0**(-0.4*(Mabs-mstar))

        mm   = np.linspace(mstar-10.,Mabs,1000)

        ll   = 10.0**(-0.4*(mm-mstar))

        lf   = phis * ll**(alpha+1) * np.exp(-ll)

        nbar = np.log10( np.trapz(lf,mm) ) - 3*np.log10(hub)

        print("%12.2f %10.2f %10.3f %10.2f"%(maglim,Mabs,L,nbar))
    

if __name__=="__main__":
    print_nbar()
    
