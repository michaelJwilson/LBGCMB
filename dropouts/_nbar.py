def get_dropoutnbar(band='g', DES=False, year=5, printit=False):
  if band in ['g', 'r', 'i', 'z']:
    if DES and band == 'g':
      '''                                                                                                                                        
          --- g-Dropouts ---                                                                                                                              
      DES year 1  23.376  26.667                                                                                                                       
      DES year 2  23.753  47.345                                                                                                                       
      DES year 3  23.973  74.237                                                                                                                        
      DES year 4  24.129  104.633                                                                                                                       
      DES year 5  24.250  138.418                                                                                                                      
      '''

      print("\n\nLoading z0 and nbar for DES-like imaging, given %s.\n\n" % band)

      ## DES-like Galaxies per sq. deg; SV equivalent to year 5.                                                                                         
      gDES  = {'Y1': 26.667, 'Y2': 47.345, 'Y3': 74.237, 'Y4': 104.633, 'Y5': 138.418, 'SV': 138.418}

      z0    = 4.0
      nbar  = gDES['Y' + str(year)]

      return  z0, nbar

    else:
      from  goldrush.specs import samplestats

      stats = samplestats(printit=printit)

      print("\n\nLoading z0 and nbar for HSC-like imaging, given %s.\n\n" % band)

  ## BX:  z ~ 2;  LBG: z ~ 3                                                                                                                               
  elif band in ['BX', 'LBG']:
    from  reddy  import  samplestats

    ## Retrieve galaxies per sq. deg. and z0 from Reddy.py                                                                                                 
    stats = samplestats(printit = printit)

    print("\n\nLoading z0 and nbar for Reddy-like imaging, given %s.\n\n" % band)

  else:
    raise ValueError("Requested type is not available for dropout p(z).")

  return  stats[band]['z'], stats[band]['nbar']


if  __name__ == "__main__": 
    print("\n\nWelcome to dropout p(z).\n\n")

    for band in ['BX', 'LBG', 'g', 'r', 'i', 'z']:
      print  get_nbar(band, DES=False, printit=False)
    
    print("\n\nDone.\n\n")
