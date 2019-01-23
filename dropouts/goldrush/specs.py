import  numpy         as      np

from    utils         import  pprint 
from    params        import  get_params
from    collections   import  OrderedDict


params = get_params()

def samplestats(printit = False):
  '''                                                                                                                                                      
  Load specifications of Reddy BX and LBG samples from tabular data and create a dictionary
  containing (interloper free) gals. per sq. deg. HSC Goldrush sample stats.
  '''

  stats                = OrderedDict()

  for band in ['g', 'r', 'i', 'z', 'y']:
    stats[band]        = OrderedDict()

  ##  Representative redshift.    
  stats['g']['z']      = 3.8
  stats['r']['z']      = 4.9
  stats['i']['z']      = 5.9
  stats['z']['z']      = 7.0

  stats['areas']       = {}
  stats['areas']['W']  = {'XMM':     28.5,  'GAMA09H': 12.4,  'WIDE12H': 15.2, 'GAMA15H': 16.6,  'HECTOMAP': 4.8,  'VVDS': 5.1}   ## deg^2
  stats['areas']['D']  = {'XMM_LSS':  2.4,  'COSMOS':   6.5, 'ELAIS_N1':  3.3, 'DEEP2_3':  5.5}                                   ## deg^2
  stats['areas']['UD'] = {'SXDS':     1.1,  'COSMOS':   1.3}                                                                      ## deg^2

  ##  Total area.
  stats['Total area']  = {}
  
  for depth in ['W', 'D', 'UD']:
    stats['Total area'][depth] = sum(stats['areas'][depth].values())  ## deg^2.

  ##  Counts in the Wide field.
  stats['g']['counts'] = {'XMM': 113582, 'GAMA9': 44670, 'WIDE': 94544, 'GAMA15': 104224, 'HECTO': 30663, 'VVDS': 23677}
  stats['r']['counts'] = {'XMM':   6371, 'GAMA9':  5989, 'WIDE':  5243, 'GAMA15':   6457, 'HECTO':  1082, 'VVDS':  1500}
  stats['i']['counts'] = {'XMM':     81, 'GAMA9':    98, 'WIDE':    36, 'GAMA15':     73, 'HECTO':    11, 'VVDS':    20}
  stats['z']['counts'] = {'XMM':      7, 'GAMA9':    16, 'WIDE':     8, 'GAMA15':     14, 'HECTO':     7, 'VVDS':    11}

  ##  Counts in the Deep field.                                                                                                                                                                                                   
  stats['g']['deepcounts']  = {'XMM': 6730, 'COSMOS': 45767, 'ELAIS': 19631, 'DEEP2': 35963}                                 ## g-band dropout at z=4.                                                                         
  stats['r']['deepcounts']  = {'XMM':  711, 'COSMOS': 6282,  'ELAIS': 612,   'DEEP2': 1498}
  stats['i']['deepcounts']  = {'XMM':    6, 'COSMOS': 64,    'ELAIS': 15,    'DEEP2': 47}
  stats['z']['deepcounts']  = {'XMM':    0, 'COSMOS': 4,     'ELAIS': 1,     'DEEP2': 5}
  

  ## Counts in the Ultra Deep field.
  stats['g']['udeepcounts'] = {'SXDS': 9916, 'COSMOS': 10644}
  stats['r']['udeepcounts'] = {'SXDS': 1209, 'COSMOS':  1990}  
  stats['i']['udeepcounts'] = {'SXDS':   36, 'COSMOS':    50}
  stats['z']['udeepcounts'] = {'SXDS':    0, 'COSMOS':     0}

  ## Linear bias as a fn. of redshift and m_UV.                                                                                                            
  stats['g']['linbias']   = {'24.0': 6.66, '24.5': 5.29, '25.0': 4.54, '25.5': 4.02}
  stats['r']['linbias']   = {'24.0': 8.47, '24.5': 7.14, '25.0': 6.21, '25.5': 5.59}
  stats['i']['linbias']   = {'25.0': 8.88}

  ## Schechter fn. fits for hsc dropouts, i.e. Table 8 of https://arxiv.org/pdf/1704.06004.pdf:
  stats['g']['schechter'] = {'phi_star': 3.040e-3, 'M_star': -20.63, 'alpha': -1.57, 'chi2_v': 6.4}
  stats['r']['schechter'] = {'phi_star': 1.060e-3, 'M_star': -20.96, 'alpha': -1.60, 'chi2_v': 1.1}
  stats['i']['schechter'] = {'phi_star': 0.540e-3, 'M_star': -20.91, 'alpha': -1.87, 'chi2_v': 0.9}
  stats['z']['schechter'] = {'phi_star': 0.438e-3, 'M_star': -20.77, 'alpha': -1.97, 'chi2_v': 2.5}

  for band in ['g', 'r', 'i', 'z']:
    ## Counts available for Wide only. 
    stats[band]['total-counts']            = sum(stats[band]['counts'].values())
    stats[band]['total-deepcounts']        = sum(stats[band]['deepcounts'].values())
    stats[band]['total-udeepcounts']       = sum(stats[band]['udeepcounts'].values())

    stats[band]['nbar']                    = stats[band]['total-counts'] / stats['Total area']['W'] 
    stats[band]['deepnbar']                = stats[band]['total-deepcounts'] / stats['Total area']['W']
    stats[band]['udeepnbar']               = stats[band]['total-udeepcounts'] / stats['Total area']['W']

    stats[band]['schechter']['phi_star']  /= params['h_100']**3.                           ## Note: h conversion is necessary.

  if printit:
    pprint(stats)
                                                     
  return  stats                                                                            ## [\phi*] = [h_100/Mpc]^3 per mag;


if __name__ == "__main__":
  print "\n*** HSC Wide survey ***"
  print "*** Total area:     %.1lf deg^2;  Limiting mag.: ~26. ***\n" % 82.6
  
  specs = samplestats()

  pprint(specs)

  print("\nDone.\n")
