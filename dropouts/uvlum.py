import numpy             as np
import astropy.constants as const


def get_uvlum(ls, Ll):
  """                                                                                                                                                                             
  Given the rest wavelength and (extincted) L_lambda                                                                                                                              
  return a measure of the UV luminosity.                                                                                                                                           
  """

  ## Rest-frame monochromatic UV luminosity defined by Sec. 6.2.3 of ./lephare/docs/lephare_doc.pdf                                                                               
  Luv  = np.trapz(Ll[(ls > 2100.) & (ls < 2500.)], ls[(ls > 2100.) & (ls < 2500.)]) * (2300 ** 2.) / (400. * const.c.value * 1e10)

  return  Luv
