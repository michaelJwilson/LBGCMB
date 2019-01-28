import  numpy              as      np

from    params             import  get_params
from    utils              import  comoving_distance
from    astropy            import  constants          as const
from    cosmo              import  cosmo
from    prep_Llls          import  prep_Llls
from    pmh                import  Pmm
from    scipy.integrate    import  simps


params = get_params()

def lensing_kernel(z):
  ## Return the lensing kernel at z.                                                                                            
  chi          = comoving_distance(z)
  chi_scatter  = comoving_distance(params['zscatter'])

  return  (3./2.)*params['om_m']*(100.*params['h_100']/const.c.to('km/s').value)**2.*(1. + z)*chi*(1. - chi/chi_scatter)

def Ckk(Pk_interps, Llls, pickle=False, zmax=params['zscatter']):
 '''
 Given a matter power spectrum interpolator and an array of L values, 
 return Ckk(L) in the extended Limber approximation.
 '''

 z           = np.linspace(0.001, zmax, 2000)
 chis        = comoving_distance(z)

 ## Stores the matter power spectrum for each of (L, z) and therefore k
 ## in the extended Limber approximation.                                                                                  
 result      = np.zeros((len(z), len(Llls)))

 for i, redshift in enumerate(z):
   ks           = (Llls + 0.5) / chis[i]
   result[i,:]  = Pmm(Pk_interps, ks, redshift)

 prefactor   = (const.c.to('km/s').value / cosmo.H(z).value)*(lensing_kernel(z) / chis)**2.

 integrand   =  prefactor[:, None] * result
 integrand  /=  params['h_100']**3.           ## account for [h^-1 Mpc]^3 of Pmm.                                                                         

 result      =  simps(integrand, dx = z[1] - z[0], axis=0)

 if  pickle:
   from  pickle  import  dump

   ## New beam/noise configurations require pickle files to be removed.
   dump(result, open("pickle/kk.p", "wb"))

 ## Ckk at each L, after integrating the lensing kernel over redshift.                                                                                   
 return  result

def var_Ckk(Llls, fsky, nkk, Pk_interps=None, samplevar_lim=False):
  ckk    = Ckk(Pk_interps, Llls)

  if samplevar_lim:
    num  = ckk

  else:
    num  = ckk + nkk

  num    =  2. * num**2.
  den    = (2. * Llls + 1.) * fsky

  return  num/den

if __name__ == "__main__":
  import  pylab  as      pl


  print "\n\nWelcome to lensing.\n\n"

  NLlls, Llls, nmodes = prep_Llls()
  zs                  = np.arange(0.0, 4.0, 0.1)
  
  kernel              = lensing_kernel(zs) / (100. * params['h_100'] * cosmo.efunc(zs))

  print "\n\nDone.\n\n"
