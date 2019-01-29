import  numpy               as  np
import  astropy.constants   as  const


def magAB(vs, mag, redshift, line = 912., depth=None):
  '''                                                                                                                                                     
  Fv for an AB mag. source of const. Fv and                                                                                                                
  given magnitude; with a Lyman break at

  lo = (1. + redshift) * line
  '''    

  dm             =  mag - 0.0
  dm            /= -2.5
  
  line          *= (1. + redshift)                                ##  Redshifted
  line           = (1e10/line) * const.c.to('m/s').value          ##  Frequency domain.

  vs            /= (1. + redshift)
  norm           = 3631. * 1e-23

  Fv             = norm * 10. ** dm * np.ones_like(vs)            ## [ergs/s/Hz].
  Fv[vs > line]  = 0.0 
  
  if depth is not None:
    dm           =  depth - 0.0
    dm          /= -2.5

    mean         = norm * 10. ** dm

    Fv          += mean * np.random.uniform(0.0, 1.0, size=len(vs))

  return  vs, Fv


if  __name__ == "__main__":
  print("\n\nWelcome to AB source.\n\n")
  
  depth                    = 29.5
  dlambda                  =  0.1

  for row in np.loadtxt("dat/sample_zm.dat"):
    redshift               = row[0]
    magnitude              = row[1]
    exposure               = row[2]

    print(redshift, magnitude, exposure)

    ls                     = np.arange(dlambda, 1.e4, dlambda)                                    ## Restframe [A].                                      
    vs                     = (1e10/ls) * const.c.to('m/s').value                                  ## ls [AA]; vs [Hz]; c [m/s.]                             

    vs, Fv                 = magAB(vs, magnitude, redshift = redshift, line = 912., depth=depth)  ## Flat Fv, Ly-lim break.
                                                                                                  ## Note: redshifts frequencies.                           

    ## vs, Fv              = redshift_spectra(vs, Fv, redshift)                                   ## redshifted spectra, observed frequencies.              

    ls                     = (1e10/vs) * const.c.to('m/s').value                                  ## Redshifted wavelengths.                               
    Fl                     = vs * Fv / ls
    
    Fl                    /= 1.e-17                                                                        ## [ergs/s/cm^2/A]

    np.savetxt("../desihub/desimodel/data/spectra/magABsource/spec-magABsource_%.2lf-LyLB_z%.1lf.dat" % (magnitude, redshift), np.c_[ls, Fl], fmt="%.6lf") 
  
  print("\n\nDone.\n\n")
