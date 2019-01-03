import  extinction
import  numpy            as      np

from    scipy.interpolate import  interp1d


def calzetti(ils, Fi, Ap=0.1, Rp=4.05, name = 'NGC6090'):
  """
  Intrinsic shape of the stellar emission, F_i(lambda), is recovered using the starburst reddening curve: 
  k'(lambda) = A'(lambda) / E_s(B-V)  with F_i(lambda) and F_o(lambda) the intrinsic and observed stellar continuum flux densities. 
  The colour excess of the stellar continuum E_s(B-v) is linked to that of nebular gas emission lines by E_s(B-v) = (0.44 \pm 0.03) E(B-V). 
  """
  Fo       = np.copy(Fi)
  ls       = np.copy(ils)

  ls      /= 10.**4.        ## Angstroms to microns. 
    
  ## E(B-V).
  Gal_Es   = {'NGC6090': 0.60, 'NGC7673': 0.41, 'NGC5860': 0.65, 'ICI586': 0.50, 'Tol1924-416': 0.02, 'Mrk66': 0.00, 'ESO185-IG013': 0.16}  
 
  E        = Gal_Es[name]
  Es       = 0.44*E
  
  kp       = 2.659*(-1.857 + 1.040/ls)
  Fo[(ls > 0.63) & (ls < 2.20)] *= 10.**(-0.4*Ap*(1. + kp[(ls > 0.63) & (ls < 2.20)]/Rp))  ## 0.63 um < ls < 2.20 um.

  kp       = 2.659*(-2.156 + 1.509/ls - 0.198/ls**2. + 0.011/ls**3.)  
  Fo[(ls > 0.12) & (ls < 0.63)] *= 10.**(-0.4*Ap*(1. + kp[(ls > 0.12) & (ls < 0.63)]/Rp))  ## 0.12 um < ls < 0.63 um.

  return  Fo

def extinct(ls, EBV=0.1, HyperZ = True, type='calzetti', atmos=False, tatmos='eso', unit='aa'):
  """
  http://extinction.readthedocs.io/en/latest/                                                                                                             
  https://github.com/kbarbary/extinction/blob/19be03f7e04ce22802c52137205aa67ae7a0a8de/extinction.pyx

  ls: wavelength in angstroms; alternative unit of 'invum' (inverse microns).                                                                           
  A_V: extinction in magnitudes at characteristic V band.                                                                                               
  Ratio of total to selective extinction, Rv = A_V / E(B - V).                                                                                         
  
  -- https://arxiv.org/pdf/1209.2152.pdf                                                                                                             
  -- https://arxiv.org/pdf/astro-ph/9911459.pdf                                                                                                      
  
  Both studies used template spectra including dust attenuation following the Calzetti law.                                                              
  Finkelstein et al. inferred an attenuation at 1500A, A1500, of 1.3 magnitudes at z=4                                                                
  and A1500 < 0.25 at z=7. In contrast, for galaxies with z = 6.5, McLure et al. found A1500=0.4,                                                     
  a value above the upper limit found by Finkelstein et al. In many observational studies the dust                                                     
  attenuation is inferred from the UV continuum slope estimated from a single colour. Bouwens et al. (2011)                                          
  measured an average UV continuum slope of -3 for galaxies at z = 7. However, this value was measured to be -2.2                                      
  when more data were collected by the HST WFC3 (Bouwens et al. 2012). This illustrates how the scarcity of high                                       
  redshift data can bias the estimation of dust attenuation.                                                                                         
      
  UV continuum slope is a poor indicator of dust attenuation, as our results are extremely sensitive to the choice                                      
  of extinction curve used as the input to the attenuation of starlight calculations. For our standard choice, a Milky Way (MW)                          
  extinction curve, galaxies get bluer when they are attenuated, which, as we explain in Section 6 is due to the presence of a                          
  bump in the extinction curve.                                                                                                             
  Usually the term dust 'extinction' refers to the attenuation of the light from a point source placed behind a screen of dust.                          
  Thus, the 'extinction' is independent of the geometry of the system.                                                                                    

  With a_\lambda = L_\lambda (attenuated) /L_\lambda (unattenuated), t_eff (\lambda) = -ln(a_\lambda).                                                   
  Magnitudes: \tau_\eff = A_\lambda / (2.5 \log_10 e).                                                                                                  
  
  MW extinction: Cardelli, Clayton & Mathis (1989) with A_V = 1.0 (Amplitude) and R_V = 3.1 (Shape).                                                     
  """

  ## Av = Rv * E(B-V); Rv = 4.05
  if HyperZ == False:
    if type == 'odonnell':
      ## MW extinction O'Donnell (1994)                                                                                                                    
      Rv  = 3.1 
      ext = extinction.odonnell94(ls, Rv * EBV, Rv)                      
                                                                                                  
    elif type == 'fitzpatrick':
      ## LMC extinction; Fitzpatrick (1999)                                                                                                                 
      Rv  = 3.10
      ext = extinction.fitzpatrick99(ls, Rv * EBV, Rv)

    else:
      Rv  = 4.05
      ext = extinction.calzetti00(ls, Rv * EBV, Rv)  ## Extinction in magnitudes at each input wavelength. 

  else:
    """
    http://webast.ast.obs-mip.fr/hyperz/hyperz_manual1/node10.html                                                                                      
    The extinction curve are expressed in  k(lambda[A]) vs lambda[A].                                                                                   
    Apply a E(B-V) correction by:                                                                                                                            
    flux_attenuated = flux_intrinsic * 10^[-0.4*k(l)*Av/Rv] 
    """

    files = {'fitzpatrick':'LMC_Fitzpatrick.dat','allen':'MW_Allen.dat','seaton':'MW_seaton.dat','prevot':'SMC_prevot.dat','calzetti':'SB_calzetti.dat'}
    Rvs   = {'fitzpatrick':                 3.10,'allen':          3.10,'seaton':           3.10,'prevot':            2.72,'calzetti':             4.05}

    files['calzetti_mod'] = 'SB_calzetti_mod.dat'             ## modified Calzetti law including contribution from 2175A bump (Massaroti et al., 2001)
    Rvs['calzetti_mod']   =                 4.05

    data       = np.loadtxt('ext-laws/' + files[type])
    extinterp  = interp1d(data[:,0], EBV*data[:,1], kind='cubic', bounds_error=False, fill_value=0.0)
    ext        = extinterp(ls)

  if atmos == True:
    files      = {'eso':'extinc_eso.dat', 'ctio': 'extinc_ctio.dat'}  ## Atmospheric extinction curves
    Rv         = 3.1

    data       = np.loadtxt('ext-laws/' + files[tatmos])
    extinterp  = interp1d(data[:,0], EBV*data[:,1], kind='cubic', bounds_error=False, fill_value=0.0)
    
    ext       += extinterp(ls)

  return 10.**(-0.4 * ext)  ## Deredden with -Av; Positive extinction values decrease flux.


if __name__ == "__main__":
  print  extinction.__file__

