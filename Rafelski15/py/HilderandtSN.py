import numpy as np
import pylab as pl


@np.vectorize
def merr(m, mstar, estar=0.2, alphab=-0.25, alphaf=0.22):
  ##  Magnitude error, 
  if m < mstar:
    return  estar * 10. ** (0.4 * (alphab + 1.) * (m - mstar))

  else:
    return  estar * np.exp(10. ** (alphaf * (m - mstar))) / 2.72 

def ferr(m, mstar, estar=0.2, alphab=-0.25, alphaf=0.22):
  flux = 10. ** -((m + 48.60) / 2.5)  ##  Assumes AB mag.  [erg/s/cm2/Hz]. 
  
  sigm = merr(m, mstar, estar=estar, alphab=alphab, alphaf=alphaf)

  ##  E.g. eqn. (6) of https://arxiv.org/pdf/1509.00870.pdf
  sigf = flux * sigm * np.log(10.) / 2.5

  return  sigf

def snr(m, mstar, estar=0.2, alphab=-0.25, alphaf=0.22, lim_snr=None):
  flux = 10. ** -((m + 48.60) / 2.5)  ##  Assumes AB mag.  [erg/s/cm2/Hz].

  ##  Depth map.
  sigf = ferr(m, mstar, estar=estar, alphab=alphab, alphaf=alphaf)

  
  if lim_snr is not None:
    sigf = np.sqrt(sigf**2. + (flux / lim_snr) ** 2.)

  snr  = flux / sigf

  return  snr


if __name__ == '__main__':
  ms = np.arange(16., 29., 0.1)
    
  for mstar in [27., 26., 25.]:
    pl.semilogy(ms, snr(ms, mstar),  label=mstar)
    pl.semilogy(ms, snr(ms, mstar,  lim_snr=100.), '--')

  pl.legend(loc=3)

  pl.ylim(0.1, 1.e3)

  pl.xlabel('m')
  pl.ylabel(r'S/N')

  pl.savefig('../plots/Hildebrandt.pdf')

    
