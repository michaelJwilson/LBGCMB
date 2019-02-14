from    __future__         import  division
from    schechterfn        import  SchechterMfn 
from    params             import  get_params
from    scipy.interpolate  import  interp1d

import  os
import  numpy  as  np


params = get_params()

def load_tables(printit = False):
  import  os
  import  pandas as pd

  root      = os.environ['LBGCMB']
  root     += '/dropouts/dat/reddy/'

  ##  Size in arcmin^2.
  ##  Notes:
  ##         In Table 1, NXX are photometric candidates, not those spectroscopically confirmed.  See https://arxiv.org/pdf/0706.4091.pdf 
  ##         In Table 3, R is the lower limit on the mag bin (typically, 0.5 in width excepth for the first which is 1.0 in width.)
  ##         Mixture of columns from Table 2 and 3, which share common row definitions.

  tab_one   = pd.read_csv(root + 'tabone.dat',   comment='#', header=None, sep=r'\s+', names=['Field', 'Size', 'NBX', 'NLBG'], engine='python')  
  tab_three = pd.read_csv(root + 'tabthree.dat', comment='#', header=None, sep=r'\s+', names=['R', 'NBX_phot', 'NBX_spec', 'NBX_int', 'NBX_fAGN',\
                                                                                              'NBX_fint', 'NLBG_phot', 'NLBG_spec', 'NLBG_int',\
                                                                                              'NLBG_fAGN', 'NLBG_fint'], engine='python')
  if printit:
    print("\n\nReddy Table one:")
    print(tab_one)
    
    print("\n\nReddy Table three:")
    print(tab_three)

  return tab_one, tab_three

def samplestats(mag=23., printit=False, h70=False):
  '''
  Load specifications of Reddy BX and LBG samples from tabular data and create a dictionary containing (interloper free) gals. per sq. deg. 
  Note:  derived from https://arxiv.org/pdf/0706.4091.pdf 
  '''

  tab_one, tab_three = load_tables(printit=printit)  
  stats              = {'BX': {'z': 2.3, 'frac_AGN': 0.03, 'schechter': {}}, 'LBG': {'z': 3.05, 'frac_AGN': 0.02, 'schechter': {}}}
  
  for survey in ['BX', 'LBG']:      
    if mag < 19.0:
      raise  ValueError('Brightest characteristics available for Reddy are R > 19.0')

    elif mag <= 22.:  
      ##  First bin is 19.0 < 22.0
      index                              = 0
      stats[survey]['Rlim']              = 22.0

    elif mag <= 25.9:
      ##  Upper limit to the bin is tab_three['R'] + 0.5;  Index where upper limit on bin is closest to provided depth. 
      index                              = np.where(np.abs(np.array(tab_three['R']) + 0.5 - mag).min() == np.abs(np.array(tab_three['R']) + 0.5 - mag))[0][0]
      stats[survey]['Rlim']              = np.array(tab_three['R'])[index] + 0.5

    else:
      raise  ValueError('Faintest characteristics available for Reddy are R < 26.0')
    
    stats[survey]['N']                   = np.array(tab_three['N' + survey + '_phot'])[:(1 + index)].sum()             
    stats[survey]['TotalArea [deg^2]']   = tab_one['Size'][tab_one['N' + survey].notnull()].sum() / 60.**2.

    ##  Note:  spectroscopic followup is not a fair sample (brighter galaxies targeted -- TBC.)  Must multiply by app. magnitude. 
    stats[survey]['_frac_interloper']    = np.array(tab_three['N' + survey + '_int']) / np.array(tab_three['N' + survey + '_spec'])

    ##  [deg^2]
    stats[survey]['nbar']                = stats[survey]['N'] / stats[survey]['TotalArea [deg^2]']

    stats[survey]['nbar_noint']          = (1. - stats[survey]['_frac_interloper']) * np.array(tab_three['N' + survey + '_phot'])
    stats[survey]['nbar_noint']          = stats[survey]['nbar_noint'][:(1 + index)].sum() / stats[survey]['TotalArea [deg^2]']
    
  ##  Schechter fn. parameterisation of the z~3 Reddy luminosity fn., Table 7 of https://arxiv.org/pdf/0706.4091.pdf 
  stats['LBG']['schechter']['phi_star']  = 1.66e-3    ## [\phi*] = [h_70/Mpc]^3 per mag for M*_AB(1700 \AA).
  stats['LBG']['schechter']['M_star']    =  -20.84
  stats['LBG']['schechter']['alpha']     =   -1.57

  ##  Schechter fn. parameterisation of the z~2 Reddy luminosity fn., Table 7 of https://arxiv.org/pdf/0706.4091.pdf                                                                                                                   
  '''
  stats['BX']['schechter']['phi_star']     = 1.74e-3     ## [\phi*] = [h_70/Mpc]^3 per mag for M*_AB(1700 \AA).                                         
  stats['BX']['schechter']['M_star']       =  -20.97
  stats['BX']['schechter']['alpha']        =   -1.84
  '''

  stats['BX']['schechter']['phi_star']     = 3.31e-3     ## [\phi*] = [h_70/Mpc]^3 per mag for M*_AB(1700 \AA).                                    
  stats['BX']['schechter']['M_star']       =  -20.60
  stats['BX']['schechter']['alpha']        =   -1.60


  if not h70:
    ##  Convert phi_star to (h_100 / Mpc)^3 per mag.
    stats['LBG']['schechter']['phi_star'] *= (10. / 7.) ** 3.
    stats['BX']['schechter']['phi_star']  *= (10. / 7.) ** 3.

    ##  Reddy (2008) Schechter fn. is in unit of h70.  MAB (1700A) - 5log10(h70) = MAB (1700A) - 5log10(h100) - 0.77                                 
    ##  Schechter fn. of M depends only on LL = M-M*, implies M* + 0.77
    stats['LBG']['schechter']['M_star']   += 0.77
    stats['BX']['schechter']['M_star']    += 0.77

    stats['LBG']['schechter']['M_star']   += 5. * np.log10(params['h_100'])
    stats['BX']['schechter']['M_star']    += 5. * np.log10(params['h_100'])

  for survey in ['BX', 'LBG']:
    if printit:
      print("\n\n%s Survey Specifications (at depth R=%.3lf):\n" % (survey, mag))

      for k, v in zip(stats[survey].keys(), stats[survey].values()):
        print("%s \t\t %s" % (k.ljust(20), v))
  
  return  stats

def Reddy_colourcut(zs, luv, u, g, r, i, z, y, type='BX', nocolourcut = False):
    if type == 'BM':
        ##  1.5 < z < 2.0                                                                                                                                   
        crit  = (g - r) >= -0.2
        crit &= (u - g) >= (g - r)     - 0.1
        crit &= (g - r) <= 0.2*(u - g) + 0.4
        crit &= (u - g) <= (g - r)     + 0.2

    elif type == 'BX':
        ##  2.0 < z < 2.5                                                                                                                                   
        crit  = (g - r) >= -0.2
        crit &= (u - g) >= (g - r)     + 0.2
        crit &= (g - r) <= 0.2*(u - g) + 0.4
        crit &= (u - g) <= (g - r)     + 1.0

    else:
        raise ValueError("Specified Reddy colour selection is not available.")

def get_pEBV(printit=False):
  names       = ['loEBV', 'hiEBV', 'p(z=2)', 'p(z=3)']

  data        = pd.read_csv("dat/reddy_tabfive.dat", sep='\s+', skiprows=0, names=names)
  data['EBV'] = 0.5 * (data['loEBV'].values + data['hiEBV'].values)

  if printit:
    for i, x in enumerate(data['EBV'].values):                                                                                                            
      print(x, data['p(z=3)'].values[i])                                                                                                                     

  return data

def get_pz(interp = True):
  zs, ns    =  np.loadtxt(os.environ['LBGCMB'] + '/dropouts/dat/reddy/pz.dat', unpack=True)
  
  ngal      =  ns.sum()
  
  dz        =  zs[1] - zs[0]

  ##  Fraction of galaxies in that bin.
  ps        =  ns / ngal
  
  ##  Probability density. 
  ps       /=  dz
  
  if interp:
    return interp1d(zs, ps, kind='linear', copy=True, bounds_error=False, fill_value=0.0, assume_sorted=False)
  
  else:
    return zs, ps
  

if __name__ == "__main__":
  import  pylab  as     pl
  from    utils  import pprint
  

  print('\n\nWelcome to the Reddy (Schechter) calculator.\n\n')

  
  Ms    = np.arange(-23., -18., 0.1)
  stats = samplestats(mag=23.0, printit=False, h70=True)

  pprint(stats)

  ##  Fig. 12 of https://arxiv.org/pdf/0706.4091.pdf
  for survey in stats:
    Phis = SchechterMfn(Ms, stats[survey]['schechter']['phi_star'], stats[survey]['schechter']['M_star'], stats[survey]['schechter']['alpha'])
    pl.semilogy(Ms, Phis, label=survey)

  pl.xlabel(r'$M_{AB}(1700\AA) - 5\rm{log}_{10}(h_{70})$')
  pl.ylabel(r'$N/mag/h_{70}^{-3}\rm{Mpc}^3$')

  pl.legend()
  pl.show()
  
  zs, ps = get_pz(interp=False)
  pl.plot(zs, ps)
  
  pl.show()
  
  print('\n\nDone.\n\n')

