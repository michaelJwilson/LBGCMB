import  numpy               as  np
import  astropy.constants   as  const


def magAB(vs, mag):
  '''                                                                                                                                                      
  Fv for an AB mag. source of const. Fv and given magnitude;
  '''    

  dm             =  mag - 0.0
  dm            /= -2.5
  
  norm           = 3631. * 1e-23                                      ##  [ergs/s/cm2/Hz].
  Fv             = norm  * 10. ** dm * np.ones_like(vs)               ##  [ergs/s/Hz].
  
  return  vs, Fv


if  __name__ == "__main__":
  print("\n\nWelcome to AB source.\n\n")
  
  magnitude                =  26.5
  dlambda                  =   0.1

  ls                       =  np.arange(dlambda, 1.e4, dlambda)       ##  Restframe [A].                                      
  vs                       = (1.e10 / ls) * const.c.to('m/s').value   ##  ls [AA]; vs [Hz]; c [m/s.]                             
 
  vs, Fv                   =  magAB(vs, magnitude)                    ##  Flat Fv
                       
  Fl                       =  vs * Fv / ls
  Fl                      /=  1.e-17                                  ##  [1.e-17  ergs/s/cm^2/A]
  
  print("\n\nDone.\n\n")
