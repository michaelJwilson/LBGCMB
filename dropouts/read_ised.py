import  os
import  array
import  numpy              as      np
import  astropy.constants  as      const

from    astropy            import  units             as u
from    redshift_spectra   import  redshift_spectra


##  Derived from read_ised of
##  https://github.com/cmancone/easyGalaxy/blob/master/ezgal/utils.py

def _read_binary(fhandle, type='i', number=1, swap=False):
  '''                                                                                                                                                     
  Reads 'number' of binary characters of type 'type' from file handle 'fhandle'                                                                            
  returns the value (for one character read) or a numpy array. 
  Set swap = True to byte swap the array after reading                                                                                                    
  '''
 
  arr = array.array(type)

  arr.fromfile(fhandle, number)

  if swap:
    arr.byteswap()

  if len(arr) == 1:
    return arr[0]

  else:
    return np.asarray(arr)

def read_ised(file = None, AgeMyr = None, printit = False, plotit = False):
  '''                                                                                                                                                
  Read a Bruzual and Charlot binary ised file.                                                                                                      
  -- input:   file (string): The name of the ised file; printit = False for verbose.                                                                
  -- returns: a tuple containing model data: (seds [ergs/s/Hz], ages [Mys], vs [Hz])                                                                
  '''
  
  if file is None:
    ##  Default to EzGal template.
    file = './SED/GAL/EZGAL/models/bc03_exp_1.0_z_0.02_chab.model'

  if not(os.path.isfile(file)):
    raise  ValueError('The model file was not found.')

  print("Reading .ised from %s" % file)

  fh    = open(file, 'rb')

  junk  = _read_binary(fh)
  nages = _read_binary(fh)

  # First consistency check.                                                                                                                             
  if nages < 1 or nages > 2000:
    raise  ValueError('Problem reading .ised; Unexpected data found for the number of ages.')
  
  else:
    print('\nNumber of ages: %d' % nages)

  ##  Read ages                                                                                                                                           
  ages  = np.asarray(_read_binary(fh, type='f', number=nages))
  ages /= 10.**9.                                                    ## Convert to units of Gyrs.                                                       

  ##  Read in a bunch of stuff that I'm not interested in but which I read like this to make sure I get to the right spot in the file.                     
  junk = _read_binary(fh, number=2)
  iseg = _read_binary(fh, number=1)

  if iseg > 0:
    junk = _read_binary(fh, type='f', number = 6 * iseg)

  junk = _read_binary(fh, type='f', number=3)
  junk = _read_binary(fh)
  junk = _read_binary(fh, type='f')
  junk = _read_binary(fh, type='c', number=80)
  junk = _read_binary(fh, type='f', number=4)
  junk = _read_binary(fh, type='c', number=160)
  junk = _read_binary(fh)
  junk = _read_binary(fh, number=3)

  ##  Read in the wavelength data.                                                                                                                         
  nvs  = _read_binary(fh)

  ##  Consistency check                                                                                                                                   
  if nvs < 10 or nvs > 12000:
    raise  ValueError('Problem reading ised file -- unexpected data found for the number of wavelengths!')
  
  else:
    print('\nNumber of wavelengths %d  ' % nvs)
  
  ##  Read wavelengths [Angstroms].                                                                                  
  ls   = _read_binary(fh, type='f', number=nvs)
  
  ##  Create an array for storing SED info.                                                                                                                
  seds = np.zeros((nvs, nages))

  ##  Now loop through and read in all the ages          
  for i in range(nages):
    junk  = _read_binary(fh, number=2)
    nv    = _read_binary(fh)
    
    if nv != nvs:
      raise  ValueError('Problem reading .ised;  Unexpected data found while reading seds.')

    seds[:,i] = _read_binary(fh, type='f', number=nvs)

    nx        = _read_binary(fh)
    junk      = _read_binary(fh, type='f', number=nx)

    if printit:
      print  "\n\nAge: %.3le" % ages[i]
      print  seds[:,i]

  fh.close()

  if AgeMyr is not None:
    seds = seds[:, np.argmin(np.abs(ages - AgeMyr))]            ##  Get the SED at that age [Lo/A].
    ages = ages[   np.argmin(np.abs(ages - AgeMyr))]            ##  Get only the SED with age closest to AgeMyr.                                             

    print("\n\nAge closest to %.3lf Myr:  %.3lf [Myr]\n\n" % (AgeMyr, ages))

  ##  Now convert the SEDS from Lo/A to Fv [ergs/s/Hz].                                                                                                     
  ## seds  *=  3.826e33 * ls ** 2.0 / (const.c.value * 1e10)    ##  [ergs/s/Hz]
  
  ls     = ls * u.AA                                            ##  Angstroms
  vs     = ls.to(u.Hz, equivalencies = u.spectral())            ##  Hertz

  ## Ll     = (seds * vs**2.) / const.c.value                   ##  [ergs/s/meter].                                                             
  ## Ll    *= 1e-10                                             ##  [ergs/s/AA].
  
  if plotit:
    import pylab as pl

    for ii, sed in enumerate(seds.T):
      if ii % 25 == 0:
        pl.loglog(ls, vs * sed, label = 'Age: %.2lf' % ages[ii] + ' Gyr.')

    pl.xlim(300., 3.e4)
    ## pl.ylim(1e13, 1e21)

    pl.xlabel(r'Wavelength [$\AA$]')
    ##  pl.ylabel(r'Flux density [ergs/s/Hz]')

    pl.legend(ncol=2, loc=4)
    
    pl.show()
    ## pl.savefig('plots/ised.pdf')

  ##  Age [Myr]; vs [Hz]; SEDS [ergs/s/Hz]; ls [A]; Ll [ergs/s/angstrom].   
  ##  Ordered by increasing frequency.   
  return  ages, vs, seds, ls, Ll 
  

if __name__ == "__main__":
  from  prep_filters  import  prep_filters
  

  print("\n\nWelcome to read-ised.\n\n")

  ## Age [Myr]; vs [Hz]; SEDS [ergs/s/Hz]; ls [A]; Ll [ergs/s/angstrom].
  ages, vs, Fv, ls, Ll = read_ised('GALAXEV/models/Padova1994/salpeter/bc2003_hr_m72_salp_ssp.ised', AgeMyr = None, printit = True, plotit=True)
  
  print("\n\nDone.\n\n")
