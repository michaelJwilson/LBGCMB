import numpy         as      np


def prep_classCls():
  import pandas             as      pd
 
  from   pandas             import  Series
  from   scipy.interpolate  import  interp1d
  from   scipy.interpolate  import  InterpolatedUnivariateSpline  as  interpolate
  from   cmb.DetectorNoise  import  DetectorNoise

  # Dimensionless total [l(l+1)/(2 * pi)] C_l's for l=2 to 4500, i.e. number of multipoles equal to 4499                                      
  #                                                                                                                                               
  # C_l^dd (deflection)        = l(l+1) C_l^phi-phi                                                                                            
  # C_l^gg (shear/convergence) = 1/4 (l(l+1))^2 C_l^phi-phi                                                                                       

  ext        =       0  ## Extrapolate if True (1). 
  splinetype = 'cubic'

  cols       = ['ell', 'TT', 'EE', 'TE', 'BB', 'PP', 'TP', 'EP']
  data       = pd.read_csv("dat/RunPB/RunPB02_cl_lensed.dat", header=None, names = cols, delim_whitespace=True)

  for x in cols[1:]:
    data[x] /= data['ell']*(data['ell'] + 1.)/(2.*np.pi)

  data['kk'] = Series(0.25*(data['ell']*(data['ell'] + 1.))**2.*data['PP'], index=data.index) # Add Cl_kappa-kappa to dataframe.                   

  spectra    = {'TT':0, 'EE':1, 'TE':2, 'BB':3}
  
  '''
  lensCl_interps = {x:interpolate(data['ell'].values, data[x].values, ext=ext) for x in spectra.keys()} # ext=0 for extrapolate; 1 to return 0. Modi=0
  lensCl_interps = {x:interpolate(data['ell'].values, data[x].values + N_inst(1.,1., data['ell'].values, x[0]), ext=ext) for x in spectra.keys()}
  lensCl_interps = {x:   interp1d(data['ell'].values, data[x].values, kind=splinetype, bounds_error=False, fill_value=0.0) for x in spectra.keys()}
  '''

  lensCl_interps  = {x:interpolate(data['ell'].values, data[x].values + DetectorNoise(data['ell'].values, 1., 1., x), ext=ext) for x in spectra.keys()}
                                                                                                                                      
  data = pd.read_csv("dat/RunPB/RunPB02_cl.dat", header=None, names = cols, delim_whitespace=True)         ## Unlensed.

  for x in cols[1:]:
    data[x] /= data['ell']*(data['ell'] + 1.)/(2.*np.pi)

  data['kk']        = Series(0.25*(data['ell']*(data['ell'] + 1.))**2.*data['PP'], index=data.index)        ## Add Cl_kappa-kappa to dataframe.                        
 
  nolensCl_interps  = {x:interpolate(data['ell'].values, data[x].values, ext=ext) for x in spectra.keys()}  ## ext=0 for extrapolate; 1 to return 0.0 on extrapolated range.

  ## Test interpolation routine. 
  ## nolensCl_interps = {x:interp1d(data['ell'].values,    data[x].values, kind=splinetype, bounds_error=False, fill_value=0.0) for x in spectra.keys()}

  return (lensCl_interps, nolensCl_interps)

def plot_class(x, xp):
  import pylab              as      pl

  from   utils              import  prefactor
  from   prep_camb          import  CAMB, Clxy
  from   cmb.DetectorNoise  import  DetectorNoise
  
  '''
  data = np.loadtxt("/Users/M.J.Wilson/work/CAMB/output/berk_lensedCls.dat")  ## External CAMB data, i.e not PyCamb.                                                                               
  pl.plot(data[:,0], data[:,1], label=r"External $TT$")                                                                                         
  '''

  pl.clf()
  
  mode = x + xp

  ells = np.logspace(0.0, 4., 5000)
  
  (lensCl_interps, nolensCl_interps) = prep_classCls() ## CLASS

  pl.plot(ells, prefactor(ells, 2)*Clxy(lensCl_interps,   ells, mode, 1, 1., 1.), label=r"CLASS lensed")
  pl.plot(ells, prefactor(ells, 2)*Clxy(nolensCl_interps, ells, mode, 1, 1., 1.), label=r"CLASS unlensed")
  
  ## CAMB
  cambx = CAMB()
  
  (lensCl_interps, nolensCl_interps) = cambx.get_Cls() 

  pl.plot(ells, prefactor(ells, 2)*(Clxy(lensCl_interps,   ells, mode, 1, 1., 1.) + DetectorNoise(ells, 1., 1., mode)), 'k-', label=r"CAMB lensed")
  pl.plot(ells, prefactor(ells, 2)*Clxy(nolensCl_interps,  ells, mode, 1, 1., 1.), 'k-', label=r"CAMB unlensed")
  
  '''
  ## Modi                                                                                                                                                 
  from   modi.nkk  import  CMB

  cmb = CMB()

  cmb.plot_TECls(x, xp, 1., 1.)
  '''
  pl.title(mode)
  pl.legend(ncol=1, loc=2)

  pl.xlim(1., 10000.)

  pl.ylim(-3.*10**-11., 3.*10.**-11.) ## TE

  pl.xscale('linear')
  pl.yscale('linear')

  pl.xlabel(r'$\ell$')
  pl.ylabel(r'$\ell (\ell + 1) C^{XY}(\ell)/2 \pi$')

  pl.savefig(r'plots/fig_class.pdf')
  

if __name__ == "__main__":
  print "\n\nWelcome to prep_class.py"
  
  plot_class('T', 'T')

  print "\n\nDone.\n\n"
