import numpy as np
import pylab as pl


def beast_nz(z):
  if (2.00 <= z) & (z <= 2.5):
    return 9.75e-4  ## (Mpc/h)^3.

  if (2.50 <= z) & (z <= 3.5):
    return 1.20e-4  ## (Mpc/h)^3.  

  if (3.50 <= z) & (z <= 4.25):
    return 1.00e-4  ## (Mpc/h)^3.  

  if (4.25 <= z) & (z <= 5.2):
    return  4.2e-5  ## (Mpc/h)^3.  

  else:
    raise ValueError('Not supported for nz at z=%.4lf' % z)

def beast_bz(z):
  if (2.00 <= z) & (z <= 2.5):
    return 2.50  ## (Mpc/h)^3.                                                                                                                          

  if (2.50 <= z) & (z <= 3.5):
    return 4.0  ## (Mpc/h)^3.                                                                                                                         

  if (3.50 <= z) & (z <= 4.25):
    return 3.5  ## (Mpc/h)^3.                                                                                                                         

  if (4.25 <= z) & (z <= 5.2):
    return 5.5  ## (Mpc/h)^3. 

  else:
    raise ValueError('Not supported for bz at z=%.4lf' % z)
