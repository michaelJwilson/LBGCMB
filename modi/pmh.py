import  os
import  math 
import  numpy              as      np

from    params             import  modi_hbias, drops_linb
from    growth_rate        import  growth_factor


## Bias taken from that measured of dropouts, or alternatively for 1e12 haloes at the stated redshift. 
## biastype    = 'dropouts' 
biastype       = 'linz'

if biastype   == 'dropouts':
  biasdict     =  drops_linb
  roundfactor  =  1.                ## Match provided redshift to the nearest in the dictionary by 'rounding'.  
                                    ## Redshifts in relevant dictionary result in a given rounding factor. 
else:
  biasdict     =  modi_hbias        ## Modi et al. measurements of linear bias for 10^12 haloes from Dark Sky.
  roundfactor  =  0.5

def rround(z):
  ##  Dropouts bias available for unit steps in redshift, modi et al. halo bias also available at half steps.
  ##  Given z, round to nearest available redshift. 
  return round(roundfactor * z)/roundfactor

def print_biasz():
  print("\nFor bias type: %s \n" % biastype)

  for z in np.arange(2.0, 7.0, 1.0):
    bestfitz   = rround(z)
    b1         = biasdict[bestfitz]['b1']

    print "At z = %.2lf, the linear bias is (1. + %.3lf)" % (bestfitz, b1)

def init_PmhPhh(cambx, fname = os.environ['LBGCMB'] + "dat/ps00_hh_DarkSky_46_z000.dat", zscatter=1100.):
  ## d: b_1, d^2: b_2, A+W: 1-loop.                                                                                                                        
  import  pandas             as      pd
  from    scipy.interpolate  import  interp1d


  splinetype              = 'cubic'

  terms                   = ['k','PZ','PA','PW','Pd','Pdd','Pd2','Pd2d2','Pdd2','Ps2','Pds2','Pd2s2','Ps2s2','PD2d','PdD2d']
  data                    = pd.read_csv(fname, header=None, names=terms, delim_whitespace=True)

  Pk_interps              = {x:interp1d(data['k'].values, data[x].values, kind=splinetype, bounds_error=False, fill_value=0.0) for x in terms}

  Pk_interps['camb_lin']  = cambx.get_lpower_interpolatekz()     ## Interpolates in both z and k.
  Pk_interps['camb_nlin'] = cambx.get_nlpower_interpolatekz()

  return  Pk_interps

def Pmm(Pk_interps, k, z, type = 'nlinear'):
  ##  Required to take an array of k values, and float argument for z.                                                                                   
  a          = 1./(1. + z)

  if (type  == 'lpt') & (z >= 1.) & (z <= 3.5):
    ## print("\n\nAssuming LPT parameters of Modi++ for Pmm prediction within (1. < z < 3.5)")

    bestfitz = rround(z)
    alpha    = modi_hbias[bestfitz]['amm']

    Ploop    = Pk_interps['PA'](k) + Pk_interps['PW'](k)

    result   = (1 - alpha*k*k/2.)*Pk_interps['PZ'](k) + Ploop

    ## result[kmax > 3.0]                      ## Impose kmax cut on LPT predictions.                                                                    
    result  *= growth_factor(a)**2.            ## Growth factor is exact in redshift.                                                                     

  elif type == 'linear':  
    ## print("\n\nAssuming CAMB linear matter power spectrum.")
    result   =  Pk_interps['camb_lin'](z, k)

  elif type == 'nlinear':
    ## print("\n\nAssuming CAMB non-linear matter power spectrum.")
    result   = Pk_interps['camb_nlin'](z, k)

  else:
    raise ValueError('Type %s not available for Pmm.' % type)

  return  result

def Pmh(Pk_interps, k, z, scross=0.):
  ##  Required to take an array of k values and float argument for z.                                                                                    
  a            = 1. / (1. + z)                                                                                                                    
  bestfitz     = rround(z)
  
  if biastype == 'Modi':
    if (bestfitz >= 1.) & (bestfitz <= 3.5):                                                                                                                
      alpha    = biasdict[bestfitz]['ahm']                                                                                                                 
      b1       = biasdict[bestfitz]['b1']                                                                                                                    
      b2       = biasdict[bestfitz]['b2']                                                                                                                    
      bs2      = biasdict[bestfitz]['bs']                                                                                                                    
      bD2      = 0.0                                                                                                                                        

      Ploop    = Pk_interps['PA'](k) + Pk_interps['PW'](k)    

      ## Read in 
      result   = (1-alpha*k*k/2.)*Pk_interps['PZ'](k) + Ploop + (b1/2.)*Pk_interps['Pdd'](k) + (b2/2.)*Pk_interps['Pd2d2'](k) + (bs2/2.)*Pk_interps['Ps2'](k)
      result  += bD2*Pk_interps['PD2d'](k) + scross

      ## result[kmax > 3.0] = 0.0               ## Impose kmax cut on LPT predictions.                                                                       

      ## NOTE: Cannot simply scale by D+^2.
      result  *= growth_factor(a)**2.           ## Growth factor is exact in redshift.                                                                     

    else:
      result   = Pk_interps['camb_nlin'](z, k)  ## For z < 1. and z > 3.5, replace with CAMB SPT results. 

  elif biastype == 'dropouts':
    b1          = biasdict[bestfitz]['b1']
    result      = (1. + b1) * Pk_interps['camb_nlin'](z, k)

  elif biastype == 'linz':
    result      = (1. +  z) * Pk_interps['camb_nlin'](z, k)

  else:
    raise  ValueError("Type is not available for Pmm.")

  return  result    

def Phh(Pk_interps, k, z, scross=0.0):
  ##  Required to take an array of k values, and float argument for z. 
  a         = 1./(1. + z)  
  
  ## Round provided z to the nearest for which a bias value is available. 
  bestfitz  = rround(z)

  if biastype == 'Modi':
    if (bestfitz >= 1.) & (bestfitz <= 3.5):
      alpha    = biasdict[bestfitz]['ahm']
      
      b1       = biasdict[bestfitz]['b1']
      b2       = biasdict[bestfitz]['b2']
      bs2      = biasdict[bestfitz]['bs']

      bD2      = 0.

      Ploop    = Pk_interps['PA'](k) + Pk_interps['PW'](k)

      ## b2*b2*Pk_interps['Pd2d2'](k): incorrect normalisation, needs divided by four on input. 
      result   = (1-alpha*k*k/2.)*Pk_interps['PZ'](k) + Ploop + b1*Pk_interps['Pd'](k) + b1**2.*Pk_interps['Pdd'](k) + b2*Pk_interps['Pd2'](k) 
      result  += b1*b2*Pk_interps['Pdd2'](k)   + b2*b2*Pk_interps['Pd2d2'](k) + bs2*Pk_interps['Ps2'](k) + b1*bs2*Pk_interps['Pds2'](k)
      result  += b2*bs2*Pk_interps['Pd2s2'](k) + bs2*bs2*Pk_interps['Ps2s2'](k) + 2*bD2*Pk_interps['PD2d'](k) + 2*b1*bD2*Pk_interps['PdD2d'](k)

      ## result[kmax > 3.0]                     ## Impose kmax cut on LPT predictions.
    
      ## NOTE: Cannot simply scale by D+^2.  
      result  *= growth_factor(a)**2.           ## Growth factor is exact in redshift. 

    else:
      result   = Pk_interps['camb_nlin'](z, k)  ## For z < 1. and z > 3.5, replace with CAMB SPT results. 

  elif biastype == 'dropouts':
    b1          = biasdict[bestfitz]['b1']       ## Stored Lagrangian bias, Eulerian is (1. + b1).
    result      = (1. + b1) * (1. + b1) * Pk_interps['camb_nlin'](z, k)

  elif biastype == 'linz':
    result      = (1. +  z) * (1. +  z) * Pk_interps['camb_nlin'](z, k)

  else:
    raise  ValueError("Type is not available for Pmm.")

  return  result


if __name__ == "__main__":
  from  prep_camb  import  CAMB


  print "\n\nWelcome."

  cambx        = CAMB()                                  

  Pk_interps   = init_PmhPhh(cambx, fname="dat/ps00_hh_DarkSky_46_z000.dat", zscatter=1100.)

  print_biasz()

  print "\n\nDone."
