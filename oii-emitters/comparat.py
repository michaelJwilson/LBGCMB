import numpy             as     np
import pylab             as     pl
import matplotlib.pyplot as     plt

from   schechterfn       import SchechterLfn
from   params            import get_params


params = get_params()

##  -- Table 7 of https://arxiv.org/pdf/1605.02875.pdf --
def _LStar(z, LStar0=41.1, betaL = 2.33):  
  result  = 10. ** LStar0
  result *= (1. + z) ** betaL

  return  result  ## [ergs/s]
  
def _PhiStar(z, PhiStar0=-2.4, betaPhi = -0.73):
  result  =      10. ** PhiStar0
  result *= (1. + z) ** betaPhi

  return  result  ## [Mpc^-3]

def _alpha(z, alpha0=-1.46):
  return  10. ** alpha0  

def oii_nbar(z, logLmin=None, printit=False):
    alpha     =            _alpha(z)   
    L_star    =            _LStar(z)  ##  [ergs/s]  
    phi_star  =          _PhiStar(z)  ##  [(Mpc)^-3] 

    phi_star /=  params['h_100']**3.  ##  [(Mpc/h)^{-3}].

    if logLmin == None:
      ##  Set lower limit on luminosity integral of L*.                                                                                                  
      Lmin   =  L_star

    else:
      Lmin   =  10. ** logLmin

    wmin     =  np.log10(Lmin / L_star)

    dw       =  0.05
    ws       =  np.arange(wmin, wmin + 10., dw)  ##  w = log10(L / L*)                                                                                       
    xs       =  np.log(10) * ws                  ##  x =    ln(L / L*)                                                                                       
    ys       =  10. ** ws                        ##  y =       L / L*                                                                                        

    nbar     =  phi_star * np.exp(-ys) * 10.**((1. + alpha) * ws)
    nbar     =  np.log(10.) * np.sum(nbar) * dw

    if printit:
      print('%.3lf  %.3lf  %.3lf  %.3le  %.3le  %.3le' % (z, alpha, L_star, phi_star, nbar))

    return  nbar                                ##  [(Mpc / h)^-3]  

if __name__ == '__main__':
  Ls = np.logspace(39., 44., 100)
  '''
  for z in np.arange(0.1, 2.4, 0.5): 
    result = []

    for L in Ls:
      result.append(SchechterLfn(L, PhiStar(z), LStar(z), alpha(z)))

    result  = np.array(result)
    ##  result /= (Ls * np.log(10))
  
    pl.loglog(Ls, result, label=z)

  ##  pl.ylim(1.e-94, 1.e-83)

  pl.legend(loc=3)

  pl.xlabel(r'$L$ [ergs/s]')
  pl.ylabel(r'$\Phi \ [\rm{Mpc}^{-3}$ per log$_{10}(L)$]')

  plt.tight_layout()
  
  pl.show()
  '''
  
  print(oii_nbar(1.0, logLmin=41., printit=False))
  
  print('\n\nDone.\n\n')
