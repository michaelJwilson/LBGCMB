from __future__ import division


def load_tables(printit = False):
  import  os
  import  pandas as pd


  root      = os.environ['LBGCMB']
  root     += '/dropouts/dat/reddy/'

  ## Size in arcmin^2.
  tab_one   = pd.read_table(root + 'tabone.dat',   skiprows=1, sep=r"\s*", names=['Field', 'Size', 'NBX', 'NLBG'], engine='python')  
  tab_three = pd.read_table(root + 'tabthree.dat', skiprows=1, sep=r"\s*", names=['R', 'NBX_phot', 'NBX_spec', 'NBX_int', 'NBX_fAGN',\
                                                                                  'NBX_fint', 'NLBG_phot', 'NLBG_spec', 'NLBG_int',\
                                                                                  'NLBG_fAGN', 'NLBG_fint'], engine='python')
  if printit:
    print "\n\nReddy Table one:"
    print tab_one
    
    print "\n\nReddy Table three:"
    print tab_three

  return tab_one, tab_three

def samplestats(printit = False):
  """
  Load specifications of Reddy BX and LBG samples from tabular data and create a dictionary containing (interloper free) gals. per sq. deg. 
  """

  tab_one, tab_three   = load_tables(printit)
  
  stats                = {'BX': {'z': 2.0, 'frac_AGN': 0.03}, 'LBG': {'z': 3.0, 'frac_AGN': 0.02}}

  for survey in ['BX', 'LBG']:      
      stats[survey]['N']                   = tab_one['N' + survey].sum() 
      
      stats[survey]['TotalArea [deg^2]']   = tab_one['Size'][tab_one['N' + survey].notnull()].sum()/60.**2.
      stats[survey]['frac_interloper']     = tab_three['N' + survey + '_int'].sum() / tab_three['N' + survey + '_spec'].sum()   

      ## [deg^2]
      stats[survey]['nbar']                =                                         stats[survey]['N']/stats[survey]['TotalArea [deg^2]']

      ## Are AGN already counted in interloper fraction?
      stats[survey]['nbar_nointerlopers']  = (1. - stats[survey]['frac_interloper'])*stats[survey]['N']/stats[survey]['TotalArea [deg^2]']

      if printit:
        print "\n\n%s Survey Specifications:\n" % survey

        for k, v in stats[survey].iteritems():
          print "%s \t\t %s" % (k.ljust(20), v)

  """
  ## Schechter fn. parameterisation of the Reddy luminosity fn., phi(L) at z = 3. 
  stats['LBG']['schechter']['phi_star']  = 1.66e-3    ## [\phi*] = [h_70/Mpc]^3 per mag for M*_AB(1700 \AA).
  stats['LBG']['schechter']['M_star']    =  -20.84
  stats['LBG']['schechter']['alpha']     =   -1.57

  ## Convert phi_star to (h_100 / Mpc)^3 per mag.
  stats['LBG']['schechter']['phi_star'] *= (7. / 10.) ** 3.
  """
  return  stats


if __name__ == "__main__":
  print '\n\nWelcome to Reddy.py'

  stats = samplestats(printit = True)

  print '\n\nDone.\n\n'

