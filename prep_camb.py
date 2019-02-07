import  camb
import  numpy              as      np

from    params             import  get_params
from    cosmo              import  cosmo
from    camb               import  model, initialpower
from    scipy.interpolate  import  interp1d


splinetype      = 'cubic'
params          =  get_params()

class CAMB:
  def __init__(self, halofit_ver = 'takahashi'):
    camb.set_halofit_version(halofit_ver)

    self.pars = camb.CAMBparams()   ## Set up a new set of parameters for CAMB                                                                              

    self.pars.set_cosmology(H0=100.*params['h_100'], ombh2 = params['om_b']*params['h_100']**2., omch2 = params['om_cm']*params['h_100']**2.)
    self.pars.InitPower.set_params(ns = params['ns'], As = params['As'])

    ## Interpolate in k, given z. 
    self.calc_linear(0.0)
    self.calc_nonlinear(0.0)
    
    self.get_background() 
    
    ## Interpolate in both k and z. 
    self.get_lpower_interpolatekz()
    self.get_nlpower_interpolatekz()
    
  def calc_linear(self, z):
    self.pars.set_matter_power(redshifts = [z], kmax=4.0)  ## Note: non-linear corrections couples to smaller scales than you want

    self.pars.NonLinear          = model.NonLinear_none    ## Linear spectrum                                                                             
    self.pkresults               = camb.get_results(self.pars)
    self.kh, self.z, self.linpk  = self.pkresults.get_matter_power_spectrum(minkh=1e-2, maxkh=2., npoints = 200)
    self.lpower_interpolatek     = interp1d(self.kh, self.linpk, kind=splinetype, bounds_error=False, fill_value=0.0)
    self.sig8                    = self.pkresults.get_sigma8()

    return  self.lpower_interpolatek

  def calc_nonlinear(self, z):
    self.pars.set_matter_power(redshifts = [z], kmax=4.0)  # Note: non-linear corrections couples to smaller scales than you want

    self.pars.NonLinear          = model.NonLinear_both

    self.pkresults.calc_power_spectra(self.pars)
    
    self.kh, self.z, self.halopk = self.pkresults.get_matter_power_spectrum(minkh=1e-2, maxkh=2., npoints = 200)
    self.nlpower_interpolatek    = interp1d(self.kh, self.halopk, kind=splinetype, bounds_error=False, fill_value=0.0)

    return  self.nlpower_interpolatek

  def get_lpower_interpolatekz(self, zscatter=params['zscatter']):
    ##  Previously kmax=10., zmax=1.1 * zscatter
    pkinterp = camb.get_matter_power_interpolator(self.pars, nonlinear=False, hubble_units=True, k_hunit=True, log_interp=True, kmax=5., zmax=1.01*zscatter)
    
    self.lpower_interpolatekz = np.vectorize(pkinterp.P) 
    
    return  self.lpower_interpolatekz  ## To be called as pkinterp.P(z, k);

  def get_nlpower_interpolatekz(self, zscatter=params['zscatter']):
    pkinterp = camb.get_matter_power_interpolator(self.pars, nonlinear=True, hubble_units=True, k_hunit=True, log_interp=True, kmax=5., zmax=1.01*zscatter)
    
    self.nlpower_interpolatekz = np.vectorize(pkinterp.P)

    return  self.nlpower_interpolatekz  ## To be called as pkinterp.P(z, k);

  def get_Cls(self):
    self.pars.set_for_lmax(10000, lens_potential_accuracy=4)

    self.clresults = camb.get_results(self.pars)
    powers         = self.clresults.get_cmb_power_spectra(self.pars, raw_cl='True') ## CMB_unit='muK': scale results from dimensionless; 
                                                                                    ## raw removes 0.5*l(l+1).  
                                                                                    ## for name in powers: print name                                     
                                                           
    lensCl         = powers['lensed_scalar']                                        ## or Total                                                             
    nolensCl       = powers['unlensed_scalar']

    ## Note L = {0,1} entries will be zero by default.                                                                                                      
    self.ell       = np.arange(lensCl.shape[0])           # Corresponding ell for CAMB                                                                        

    ## "Hu 2002: under the assumption of parity invariance:          Cl_TB = Cl_EB = 0.0.                                                                  
    ##  While, in the absence of grav. waves, lensing and vorticity: Cl_BB = 0.0"                                                                          
    self.spectra    = {'TT':0, 'EE':1, 'BB':2, 'TE':3}
    
    self.lensCl_interps   = {x:interp1d(self.ell,  lensCl[:,self.spectra[x]],kind=splinetype,bounds_error=False,fill_value=0.0) for x in self.spectra.keys()}
    self.nolensCl_interps = {x:interp1d(self.ell,nolensCl[:,self.spectra[x]],kind=splinetype,bounds_error=False,fill_value=0.0) for x in self.spectra.keys()}

    # Cls for lensing potential.                                                                                                                           
    ppCl               = self.clresults.get_lens_potential_cls(self.ell.max(), raw_cl='True')[:,0]
    kkCl               = 0.25*(self.ell*(self.ell + 1.))**2.*ppCl
    kkCl_interp        = interp1d(self.ell, kkCl, kind=splinetype, bounds_error=False, fill_value=0.0)

    # Add kappa.
    self.spectra['kk']        = 4                                                                                                                          
    self.lensCl_interps['kk'] = self.nolensCl_interps['kk'] = kkCl_interp

    return (self.lensCl_interps, self.nolensCl_interps)

  def get_background(self, z = 0.0):
    self.backresults = camb.get_background(self.pars)
    self.chistar     = self.backresults.conformal_time(0) - model.tau_maxvis.value
    self.zstar       = self.backresults.get_derived_params()['zstar']

    ## print 'Derived parameter dictionary: \n'
    ## print self.backresults.get_derived_params()
    
    ## Derived parameter dictionary:
    ## {'rdrag': 146.97861407805524, 'rstar': 144.1990827302238,       'age': 13.740500129586403, 
    ##  'kd': 0.14069704144187056,   'zdrag': 1059.2079162597656,      'thetaeq': 0.8080688332825077, 
    ##  'zstar': 1090.5690167257048, 'thetarseq': 0.44700138490487534, 'zeq': 3441.168556293452, 
    ##  'thetad': 0.161748306743081, 'thetastar': 1.0445694817769533,  'keq': 0.010502751815052997, 
    ##  'DAstar': 13.804642510225527
    ## }
    
  def get_dists(self, z):
    return (backresults.comoving_distance(z), backresults.angular_diameter_distance(z), backresults.luminosity_distance(z))

  def chi2z(self, chis, nz_step=150, zmax=10000):
    return self.backresults.redshift_at_comoving_radial_distance(chis, nz_step=150, zmax=10000) ## [Mpc]. 
    
  def print_pars(self):
    print(self.pars)

## Wrapper for CAMB and CLASS Cls to add detector noise. 
def Clxy(Cl_interps, ell, alpha='TT', thetab=1., DeltaT=1.): # u.arcmin                                                 
  from   cmb.DetectorNoise  import  DetectorNoise

  
  result     = Cl_interps[alpha](ell)
  
  result    += DetectorNoise(ell, thetab, DeltaT, type=alpha)         # Add detector noise.                                   
  ## result +=        N_inst(thetab, DeltaT, ell, alpha[0])           # Will not return 0.0 for TE, but "should". 
  
  return  result


if __name__ == "__main__":
  print("\n\nWelcome.")

  print("\n\nCreating a class instance of CAMB.")

  cambx = CAMB()

  (lensCl_interps, nolensCl_interps) = cambx.get_Cls()  ## CAMB; CLASS:  prep_classCls().                                                               
  
  print('\n\nDone.\n\n')
  
