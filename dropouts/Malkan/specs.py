import  os
import  numpy        as      np
import  pylab        as      pl

from    collections  import  OrderedDict
from    params       import  get_params


params  = get_params()

def samplestats(mag=25.5):
  stats                         = {'Malkan': {'z': 3.0}, 'Reddy': {}, 'Sawicki': {}, 'Bian': {}}

  ## Table 3 of Malkan et al. (2017), https://arxiv.org/pdf/1711.04787.pdf
  ## Note:  log10 |phi* / Mpc^-3|
  stats['Malkan']['schechter']  = {'phi_star': -2.73, 'M_star': -20.86, 'alpha': -1.78}  
  stats['Reddy']['schechter']   = {'phi_star': -2.77, 'M_star': -20.97, 'alpha': -1.73}
  stats['Sawicki']['schechter'] = {'phi_star': -2.77, 'M_star': -20.90, 'alpha': -1.43}
  stats['Bian']['schechter']    = {'phi_star': -2.97, 'M_star': -21.11, 'alpha': -1.94}

  for band in stats:
    stats[band]['schechter']['phi_star']  =  10. ** stats[band]['schechter']['phi_star']
    stats[band]['schechter']['phi_star'] /=  params['h_100']**3.                           ## Note: h conversion is necessary.                                                       
                                                                                           ## [\phi*] = [h_100/Mpc]^3 per mag; 
  ##  Add in the Rc band measurements.
  dat = np.loadtxt(os.environ['LBGCMB'] + '/dropouts/Malkan/dat/Table2.dat')

  stats['Malkan']['Rlim'] = mag
  stats['Malkan']['nbar'] = dat[:,3][dat[:,1] <= mag].sum() * 60. * 60.                    ## N(m < mag) per sq. deg.
  
  return  stats


if __name__ == '__main__':
    print("\n\nWelcome.\n\n")

    stats = samplestats()
    
    for key in stats:
      for skey in stats[key]:
        print key, skey, stats[key][skey]
    
    print("\n\nDone.\n\n")
